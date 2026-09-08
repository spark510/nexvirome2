# 숙주 경쟁 판정·NMF 해석·graph 재구성 구현 계획

작성: 2026-09-07. 이 문서는 구현 예정 계약이다. 아래 새 CLI와 파일은 아직 구현 완료를 뜻하지 않는다.
2026-09-08 추가: bounded PE 기반 `reconstruct run` prototype과
`workflows/reconstruction.smk`를 구현했다. 공통 anchor 경로 교체, 좌표표,
원본 contig 경쟁 정렬과 provisional 결과 재정렬을 지원한다.
여러 특징의 통합 경로 점수·read 검색 최적화·phasing은 아직 미구현이다.
실제 현재 계약은 [reconstruction.md](reconstruction.md)를 참조한다.
기존 read preparation, multiview, correction, graph, runtime 모듈을 확장한다.
외부 도구의 실제 실행 검증과 생물학적 성능 검증을 코드 테스트와 구분한다.

## 1. 목표와 범위

1. 숙주 판정표를 자동 생성하고, 명확한 숙주 fragment만 선택적으로 제거한다.
2. NMF 성분별 독립 생물학적 근거를 요약하고 불확실성을 표시한다.
3. 기존 graph에 존재하는 대체 경로를 검증해 국소 재구성 후보를 만들고 선택적으로 적용한다.
4. 여러 분기 사이를 실제 fragment 근거로 연결할 수 있을 때만 국소 haplotype을 구성한다.

독립 assembler, 사라진 graph 경로 복구, 근거 없는 gap 채우기, NMF 성분=종/strain이라는 해석,
장거리 근거 없는 완성 strain genome 출력은 이번 구현 범위에 넣지 않는다.
현재 split-only 교정은 별도 비교군으로 계속 제공한다.

## 2. 공통 설계와 코드 구조

```text
src/nexvirome2/
  host/                 reference, mapping, evidence, policy, service
  feature_integration/  schema, builders, service
  component_annotation/ evidence, aggregation, policy, service
  reconstruction/       candidates, support, conflicts, spelling, service
  phasing/              linkage, blocks, service
  domain/               위 기능의 공통 record와 설정
```

각 디렉터리를 Python package로 만들되 파일 분리는 책임과 크기에 따라 결정한다.
상태 없는 계산은 함수로 두고, dataclass로 ReferenceBundle, FragmentEvidence,
ComponentEvidence, ReconstructionCandidate, CoordinateSegment, PhaseBlock을 표현한다.
정렬 실행 backend와 판정 policy에만 Protocol을 사용해 실제 도구와 fixture 구현을 교체한다.
기존 CommandRunner/RunContext, FASTA·SAM·TSV 입출력을 재사용한다.

공통 계약:
- 좌표는 0-based half-open, strand는 별도 필드. ORF·reference·graph 좌표계를 혼용하지 않는다.
- 모든 입력과 DB·도구·정책·특징 schema의 버전/checksum을 기록한다.
- missing, measured-zero, no-hit, not-searched를 구분한다.
- read 이름의 출처나 simulation truth는 운영 판정 입력으로 받지 않는다.
- 불완전한 외부 도구 출력, 제한에 걸린 탐색, 실패한 실행을 성공/음성 근거로 바꾸지 않는다.
- 생물학적 판정 설정은 개발 panel에서 고정한 뒤 독립 평가한다.

## 3. 자동 숙주 경쟁 판정

### 3.0 FASTQ QC 선행 단계 — 2026-09-08 추가

실험 read 입력은 Illumina paired-end R1/R2 FASTQ 또는 FASTQ.gz로 고정한다.
Reference는 FASTA이며 single-end/long-read는 별도 지원 범위다.
현재 reads prepare는 형식·길이·평균 quality·N 비율 검사만 수행하고 adapter trimming은 하지 않는다.
새 reads qc backend를 숙주 경쟁 판정 앞에 추가한다. 아래는 구현 계획이며 현재 연결된 기능이 아니다.

```text
raw FASTQ → 형식/pair 검사 → fastp 또는 Cutadapt → QC FASTQ
         → host classify → 선택적 host subtraction → assembly
```

일반 shotgun 기본 backend 후보는 fastp다. Adapter 처리, quality/길이 filtering과
HTML/JSON 보고서를 제공한다. PE overlap 기반 adapter trimming과 추가 adapter 자동 검출은
구분해서 설정하고, 알려진 adapter가 있으면 library metadata에서 명시한다.
[fastp 공식 문서](https://github.com/OpenGene/fastp)

특정 primer, anchored/linked adapter 처리가 필요한 library에는 Cutadapt backend를 제공한다.
R1/R2 adapter, 최소 overlap, 허용 오류와 pair filtering 정책을 명시한다.
[Cutadapt 공식 문서](https://cutadapt.readthedocs.io/en/stable/guide.html)

초기에는 두 backend 중 하나만 사용하며, 연속 적용은 목적과 각 단계 역할이 명시된 경우만 허용한다.
Amplicon은 primer 제거를 지원하더라도 shotgun coverage/NMF 실험과 동일한 library로 취급하지 않는다.
Quality threshold는 도구별 의미가 다르므로 현재 min_mean_quality=20을 fastp의 base-quality
threshold에 그대로 대입하지 않는다. 길이 50 bp는 기존 시작값으로 두고 read 길이·assembly 조건별 검증한다.
Sliding-window trimming, poly-G/poly-X, 저복잡도 제거, overlap correction, deduplication,
read merging은 명시적으로 정책을 기록한다. 특히 merging/correction/deduplication과
저복잡도 제거는 기본 비활성화 계획이다. Poly-G는 장비/library에 따라 설정한다.

원본 pair ID를 보존하고 trimmed/retained 길이와 제외 사유를 추적한다.
한쪽 mate만 통과한 read는 assembly용 paired 출력에 섞지 않고 별도 보존한다.
원본 FASTQ를 보존하며 도구별 failed/unpaired 출력 형식은 공통 pair 감사표로 정규화한다.
출력: clean_R1/R2.fastq.gz, QC JSON/HTML 또는 도구 보고서, pair lineage, summary, manifest.
Host metadata의 checksum은 이 clean FASTQ를 기준으로 한다.
기존 reads prepare는 QC 완료 경로에서 중복 quality filtering을 하지 않도록 validation/subtraction
모드를 분리하되 기존 기본 동작은 유지한다. 필터 생략 여부와 upstream QC manifest를 기록한다.

Snakemake에 QC rule, doctor에 선택 backend 검사, Conda 환경에 fastp/Cutadapt를 추가하고
conda-forge·bioconda만으로 Linux 의존성 해결과 실제 tool smoke test를 실시한다.
현재 환경 파일에는 아직 이 두 도구를 추가하지 않았다.
ART 기본 benchmark는 adapter 오염을 가정하지 않고 QC 전후 arm을 분리한다.
Adapter/저품질 오염 fixture를 별도로 만들어 과도한 trimming, viral pair 손실,
잔여 adapter, mate 일치, 빈 출력, variable read length, gzip 처리를 검증한다.
NMF/assembly 비교군에는 동일 QC 출력을 공급해 QC 효과와 교정 효과를 분리한다.
구현 순서는 공통 schema 직후, host classify 구현 전에 둔다.

### 3.1 Reference bundle

입력: host FASTA, viral competitor FASTA, 선택적 plasmid/decoy FASTA,
각 서열의 역할·accession·taxonomy metadata, 선택적 host ERV/공유 영역 annotation.
RefSeq Virus는 바이러스 경쟁군으로 사용하며, 별도의 숙주 reference가 필요하다.
대상 숙주 species/assembly는 프로젝트에서 임의 확정하지 않고 실행 설정으로 요구한다.

서열 ID를 역할별 namespace로 재기록해 충돌을 방지하고 원본 대응표를 보존한다.
원본 FASTA는 유지한다. ERV annotation은 host의 좌표 interval로 표시하며 중복 reference
추가만으로 새로운 경쟁 근거가 생긴 것으로 세지 않는다.
sharing atlas와 host 고유 구간 정보를 제거 판정의 보조 근거로 사용한다.
식별력 평가용 mask와 보존 서열의 전체 정렬 근거를 모두 기록한다.

출력: combined FASTA, targets.tsv, shared/ERV intervals, Bowtie2 index, bundle manifest.
인덱스 캐시는 모든 source checksum, index 옵션, Bowtie2 버전으로 식별한다.

### 3.2 Mapping과 fragment 단위 판정

초기 backend는 기존 Conda 환경의 Bowtie2 + SAMtools를 사용한다.
동일 combined index와 동일 scoring 설정에서 경쟁 정렬한다. 서로 다른 설정의 점수를 비교하지 않는다.
기본은 PE end-to-end이며 local alignment는 별도 정책으로 분리한다.
primary/secondary와 AS/XS, 실제 aligned bases, mate 방향·거리, coverage를 보존한다.

대용량 실행은 두 단계로 구성한다.
1. 제한된 경쟁 hit를 수집해 후보를 분류한다. hit 수 제한에 도달하면 competition_truncated로 표시한다.
2. 제거 후보 중 경쟁 정보가 불완전한 fragment만 더 깊게 재정렬한다. 추가 탐색이 예산을 넘으면 보존한다.

MAPQ만으로 숙주를 확정하거나, 정렬기가 모든 후보를 찾았다고 가정하지 않는다.
양쪽 mate가 숙주 구간에 일관되게 정렬되고, 경쟁군과 비교해 충분한 점수 차이와
host 식별 가능한 지지 구간을 가진 경우에만 host_supported 후보로 둔다.
정렬 범위, 점수 차이, informative bases의 임계값은 개발 fixture/panel에서 결정한다.
반대쪽 mate가 viral/plasmid를 지지하거나 ERV/공유 영역만으로 설명되면 ambiguous로 남긴다.
viral DB에 hit가 없다는 사실은 host-positive 근거가 아니다.

기존 adapter와 호환되는 statuses:
- host_supported: 제거 조건 충족.
- nonhost: 바이러스/plasmid 등 경쟁군 지지. 세부 역할은 별도 필드.
- ambiguous: 충돌, 공유 영역, 불완전한 경쟁 탐색, 낮은 식별력.
- unclassified: 충분한 정렬 근거 없음.

### 3.3 연결과 산출물

계획 CLI:
```text
nexvirome2 host reference --host ... --competitors ... --metadata ... --config ... --output ...
nexvirome2 host classify --r1 ... --r2 ... --bundle ... --config ... --output ...
```

출력: host_calls.tsv, host_metadata.json, fragment_evidence.tsv, mapping BAM, summary.json.
판정표는 모든 입력 pair를 한 번씩 포함하고 metadata에 기존 read checksum 계약을 유지한다.
기존 reads prepare의 --host-calls/--host-metadata로 바로 연결한다.
제거 활성화는 별도 설정이며 초기에는 report-only가 기본이다.
raw 입력에 대한 판정표를 raw 입력의 reads prepare에 전달해 checksum이 어긋나지 않게 한다.
QC를 먼저 했다면 QC 출력 FASTQ의 checksum으로 새 판정표를 만든다.

판정표 조회와 read-ID 중복 검사는 현재 메모리 dict/set 방식에서 디스크 기반 join으로 확장한다.
paired reads는 항상 함께 보존/분리하고 제거 reads도 출력한다.

## 4. 특징 자동 통합과 NMF 성분 해석

### 4.1 특징 통합

계획 CLI: nexvirome2 features integrate --coordinates ... --inputs ... --config ... --output ...

공통 contig/window 좌표를 기준으로 다음 block을 기존 multiview long TSV로 통합한다.
- biological sample별 coverage. technical replicate는 사전에 집계한다.
- k-mer/tetranucleotide 구성. codon usage는 유효 ORF가 있는 경우에만 측정한다.
- gene/domain family, viral hallmark, host 고유 정렬, 구조 annotation.

출력: features.tsv(entity,block,feature,value), coordinates.tsv, feature_schema.json,
feature_provenance.tsv. 측정 범위·DB·검색 여부와 누락 이유를 함께 기록한다.
길이와 유전자 수 차이를 보정하고 중복 domain/HSP를 합쳐 같은 근거의 중복 투표를 막는다.
window 크기, block 정규화, missing 처리와 threshold를 설정으로 노출한다.
기존 NMF의 observed-mask와 비음수 계약을 유지한다.

### 4.2 성분별 근거 집계

계획 CLI: nexvirome2 components annotate --model ... --coordinates ... --evidence ... --config ... --output ...

NMF loading을 가중치로 사용하되 겹치는 windows나 한 contig의 많은 ORF가
성분 전체 판단을 지배하지 않도록 contig/근거 family 단위로 먼저 집계한다.
절대 지지량, 지지 contig 수, 반대 근거, 검색 가능한 범위, seed 안정성을 함께 출력한다.

라벨은 host-supported, virus-supported, plasmid-supported, mixed, unknown으로 시작한다.
conserved_signal은 별도 축으로 출력한다. 보존성은 생물학적 출처 종류가 아니다.
host+virus 근거의 공존은 통합 바이러스 등일 수 있으므로 자동으로 chimera라고 선언하지 않는다.
성분 라벨을 개별 contig의 확정 라벨로 일괄 전파하거나 host 제거 조건으로 직접 사용하지 않는다.

같은 domain 특징을 NMF 입력과 해석에 썼으면 descriptive evidence로 표시한다.
이를 독립 검증으로 세지 않고, 별도 truth 또는 학습에 쓰지 않은 근거로 평가한다.
loading과 가중 점수는 확률이 아니다. 별도 calibration 없이는 confidence probability로 표시하지 않는다.

출력: component_annotations.tsv, component_evidence.tsv, entity_annotations.tsv,
conflicts.tsv, interpretation.json. source evidence까지 추적 가능해야 한다.

### 4.3 새 데이터 적용

기존 multiview에 transform을 추가한다. 개발 데이터에서 저장한 feature 순서·scale·H를
고정하고 새 데이터의 W만 추정한다. schema 불일치와 측정 불가능 feature를 명시한다.
평가 데이터를 포함해 scale/H를 다시 학습하는 경우 transductive 실험으로 별도 보고한다.

## 5. 국소 graph 경로 재구성

### 5.1 후보와 read 접근 최적화

Graph.spell, LocalPaths, 현 pair evidence를 재사용한다.
먼저 branch/repeat 후보, node 방향, overlap, 원본 contig 구간, 대체 path를 추출한다.
국소 재구성 1차는 양 끝이 원본 경로의 동일한 anchor로 연결되는 bounded 후보를 지원한다.
다른 exit로 끝나면서 원래 suffix와 연결되지 않는 후보는 원본 contig 대체로 적용하지 않는다.
그 경우 별도의 국소 후보 서열만 출력하거나 unresolved로 남긴다.

node 수·경로 수·서열 길이·탐색 시간 제한을 두고 초과 사유를 기록한다.
원본 경로와 경쟁 경로를 함께 출력한다. alternate entrance도 경쟁군에 포함한다.

각 후보마다 전체 reads를 재정렬하는 비용을 줄이기 위해 graph segment 전체 또는
배치 경로 library에 공통 mapping을 수행하고 read-name retrieval index를 만든다.
후보 node 및 경쟁 flank에 정렬되는 fragment의 양쪽 mate를 가져온다.
reference mapping에서 누락되는 junction-spanning reads를 위해 path-junction seed 검색을 보완한다.
fixture에서 전체 read 재정렬 대비 후보 회수율을 검증하기 전에는 전체 정렬 backend를 유지한다.

### 5.2 경로 선택 정책

각 경로에 독립 PE 지지, 반복/다중 정렬 비율, sample coverage 관계,
composition 유사도, NMF 관계, 생물학적 충돌을 연결한다.
기존 paths rank를 단순 baseline으로 사용하고 missing을 음성 점수로 바꾸지 않는다.

재구성 허용에는 최소한 다음이 필요하다.
- 모든 변경 junction에 실제 식별 가능한 fragment 지지.
- 원본 및 모든 열거된 경쟁 경로 대비 충분한 우위.
- 대체 경로를 제외시킨 탐색 제한이나 정렬 누락이 판단을 왜곡하지 않음.
- graph 방향/overlap과 양쪽 anchor의 연결 일관성.

초기 min_fragments=3, alternative_fraction=0.9, score_margin=1을 개발 시작값으로 쓰되
충분한 재구성 보증이라고 간주하지 않는다. 여러 분기를 통과하는 경로는 국소 지지의 합만으로
허용하지 않고 분기 간 linkage가 필요하다. NMF/생물학적 점수만으로 PE 근거 부족을 대체하지 않는다.

결정은 retain, split, reconstruct, unresolved로 분리한다.
실패한 재구성 후보를 이유로 자동 분할하지 않고 기존 split policy를 독립 적용한다.

### 5.3 충돌 해결과 서열 출력

변경 후보를 모아 shared node, 원본 좌표 중첩, 서로 다른 suffix 선택 등의 충돌을 검사한다.
1차는 서로 독립된 후보만 적용하고 충돌 집합은 unresolved로 남긴다.
단순 점수 우선으로 하나를 강제 선택하지 않는다.

Graph.spell로 실제 graph 서열만 출력한다. 임의 consensus base/gap filling은 하지 않는다.
원본 contig를 덮어쓰지 않고 reconstructed.fasta와 유지/분할 결과를 별도 보존한다.
완전히 동일한 결과는 서열과 경로 provenance를 보존하며 중복 정리한다.

재구성은 기존 split 좌표표와 달리 새 node 서열을 포함할 수 있으므로 schema를 분리한다.
좌표표: output_id/output_start/output_end, graph_node/node_start/node_end/strand,
선택적 original_contig/original_start/original_end, operation.
기존 원본에 없는 구간은 original 좌표를 null로 두며 억지로 일대일 대응시키지 않는다.

출력: candidates.gfa/tsv, path_evidence.tsv, decisions.tsv, conflicts.tsv,
reconstructed.fasta, graph_coordinate_map.tsv, provenance.json.
재구성 FASTA로 reads를 재정렬해 각 변경 junction의 지지와 coverage를 재검사한다.
이는 입력 reads와의 일관성 확인이며 독립적인 정확도 증명은 아니다.

계획 CLI: nexvirome2 reconstruct --contigs ... --graph ... --paths ... --r1 ... --r2 ... --config ... --output ...
기본은 후보 생성/report-only이고 검증된 정책에서만 apply를 활성화한다.

## 6. 국소 phasing과 strain 분리의 경계

재구성 기능이 검증된 뒤 별도 CLI인 nexvirome2 phase로 진행한다.
fragment가 함께 지지하는 branch allele 조합으로 linkage graph를 만든다.
branch 간 linkage가 있는 연결 성분만 phase block으로 정의한다.
독립적인 fragment 수, 상충 조합과 multimapping을 평가해 복수 haplotype 후보를 보존한다.
연결되지 않은 block은 coverage/NMF 유사성만으로 하나의 strain으로 잇지 않는다.
read 길이보다 먼 구간의 phasing은 실제 장거리 근거가 없으면 unresolved다.
출력은 phase_blocks.tsv, linkage.tsv, local_haplotypes.fasta와 미해결 구간이다.
기능 이름이나 보고서에서 이를 완성 strain genome 복원으로 표현하지 않는다.

## 7. Workflow·환경·실행 기록

host.smk: bundle → competitive mapping → host calls → 기존 reads prepare.
extensions.smk: feature integration → NMF fit/transform → component annotation.
reconstruction.smk: candidates → read retrieval → evidence → decisions → spelling → remapping/evaluation.
phase는 재구성 이후 선택적으로 연결한다.

DB, 모델, host reference는 경로/checksum으로 선언하고 설치 시 자동 다운로드하지 않는다.
Linux x86_64, 기존 Bowtie2/SAMtools/BLAST/Python/Snakemake를 우선 사용한다.
Conda 채널은 conda-forge·bioconda만 사용한다. 최초 실제 서버 검증 후 버전을 고정한다.
단계별 config, revision, 명령, tool version, seed, 자원 사용량, 실패 로그를 보존한다.
캐시는 입력/정책/도구 hash를 모두 확인하고 다른 설정의 결과를 재사용하지 않는다.

## 8. 테스트와 실험 설계

단위/통합 fixture:
- host 고유, viral 고유, shared-only, ERV, 한쪽 mate만 host, viral DB에 없는 원본.
- 동점/비정상 점수, MAPQ unavailable, 정렬 제한, secondary/supplementary, mate 누락.
- NMF component 번호 permutation, missing vs zero, 중복 domain, window 중복, 상충 근거.
- frozen transform 재현성과 feature 순서 불일치, 학습 특징을 독립 evidence로 오인하는 경우.
- graph 역방향, overlap, 반복 node, 원형 junction, alternate entrance, 연결 안 되는 exit.
- overlapping rewrite, 동일 결과 중복, 새 node 좌표, 지원 없는 gap, 분리된 phase blocks.
- truth 인자 없이 실행 가능, 실패 manifest, 경로/샘플별 checksum, disk cache 정리.

실제 도구 smoke test: 소규모 FASTA/PE fixture로 Bowtie2 인덱싱부터 workflow 끝까지 실행한다.
이후 자연 RefSeq와 별도 host/plasmid reference에서 개발·평가 panel을 구성한다.
근연 strain/공유 source가 양쪽 panel로 새지 않도록 genome/유사도 cluster 단위로 분리한다.
기술 seed는 독립 표본 수로 세지 않는다. species holdout과 viral competitor DB 누락 실험을 포함한다.

단계별 ablation을 먼저 실시하고 그 다음 필요한 조합만 통합 평가한다.
- host: 제거 없음 vs 자동 제거. QC는 동일하게 고정.
- NMF 해석: 기존 독립 근거 baseline vs NMF 성분 집계. 학습에 쓴 annotation을 정답으로 재사용하지 않음.
- 재구성: 원본, split-only, PE-only 재구성, 가중 특징, NMF 추가, 각 분할 수에 맞춘 random split.
- 전체 비교: MEGAHIT 원본과 metaSPAdes 원본도 유지. 조립 조건과 resource budget을 기록.

측정:
- host 제거율, viral pair 손실률, 공유/미지/ERV 하위군별 손실, 최종 viral genome recovery.
- component precision/recall과 unknown 비율, conflict, seed/샘플 변화 안정성.
- 새로 추가한 잘못된 junction 수, 수정한 junction 수, correctly contiguous length,
  원본별 recovery, 중복 복원, strain switch, unresolved 비율.
- 시간, peak RSS, temporary disk, 총 정렬 횟수/처리 read 수.

## 9. 실행 순서와 완료 조건

| 순서 | 산출물 | 완료 기준 |
|---|---|---|
| 1 | 공통 schema·fixture·설정 | CLI 계약 및 실패/좌표/출처 테스트 |
| 2 | host bundle·판정·read adapter 연결 | 자동 calls 생성, 보존 정책, 실제 mapping smoke test |
| 3 | 특징 통합·성분 해석·transform | 동일 좌표/결측/중복 계약, 독립 근거 평가 가능 |
| 4 | graph 후보·공통 read 검색 | 경쟁군 포함, 전체 재정렬과 후보 회수율 비교 |
| 5 | 국소 재구성·좌표표·충돌 처리 | 모든 변경 junction 검증, source fixture에서 정답 복원 |
| 6 | Snakemake 통합·Linux 환경 고정 | 작은 실제 도구 end-to-end 실행 및 재현 |
| 7 | 개발 tuning·정책 freeze·독립 평가 | 사용 모듈을 성능/손실/자원 근거로 결정 |
| 8 | 국소 phasing | linkage가 있는 block만 복원, switch error 평가 |

초기 기능은 모두 연구/선택 기능으로 배포한다. 기존 freeze/select/execute의 허용 목록은
새 계약과 설정 fingerprint, 정책 검증을 추가한 뒤에만 확장한다.
host 제거 허용 손실률과 component 라벨 채택 기준은 개발 전에 설정 파일에 명시한다.
개발 시작 제안값은 host 처리의 viral pair 손실률 0.1% 이하이며, 전체 평균 외에
미지/공유/ERV 등 하위군 결과와 신뢰구간을 확인한다. 이는 보증이나 이미 달성한 결과가 아니다.
재구성의 최초 목표는 원래 계획의 확정 오류 연결 50% 감소·recovery 감소 2 percentage points 이내를
유지하면서 새로 추가하는 오류 연결을 별도 제한하는 것이다. 구체 허용치와 최소 표본 수는
개발 panel 설계 때 고정하고 heldout 결과를 본 뒤 완화하지 않는다.
baseline 오류가 없으면 개선률을 주장하지 않고 재현되지 않은 범위와 한계를 보고한다.
효과가 없거나 손실이 큰 모듈은 report-only/연구 기능으로 남기며 최종 기본 파이프라인에 넣지 않는다.
