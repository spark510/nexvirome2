# NexVirome2 통합 구현 계획 및 요구사항 대조

갱신: 2026-09-07. 이 문서는 앞으로 구현할 작업의 기준이다. 실행 가능한 현재 기능은
[implementation.md](implementation.md), 확장 기능은 [extensions.md](extensions.md), 검증 이력은
[validation.md](validation.md)에 따로 기록한다. 아래 `미구현`·`부분 구현`·`후속`은 완료로 간주하지 않는다.
추가로 개발 조건별 mask 후보 선택, 결측값을 제외한 혼합 특징 NMF와 교정 gate,
보존 read 근거 재사용, paired-read QC/host-call adapter, 근접 서열 catalog 및 선택적 coassembly를
구현했다. 구체적인 입력·제한·실행법은 [remaining_extensions.md](remaining_extensions.md)에 기록한다.
2026-09-07에 taxonomy/shared-seed atlas, mask import/propose/export/평가, annotation 좌표 연결,
joint 집계, discovery 근거 카드, exact catalog, 그룹 holdout, 생물학적 다중 샘플 설계,
graph 후보 보존, 구조 표현 batch cache, geNomad 결과 adapter 및 운영 freeze/select/execute를 추가했다.
원형 seed의 원본 좌표 분리와 BLAST 기반 근사 공유 atlas도 추가했으며, 확장 workflow는
원본·교정 5종·각각의 random control 총 11개 arm의 실제 평가 경로를 제공한다.
실행법과 세부 제한은 [masking_pipeline.md](masking_pipeline.md)를 따른다. 아래 요구사항 표의
상태는 최초 대조 시점의 기준선이며, 이 추가 구현으로 R02–R07/R10–R14/R16–R24의 일부가 진전됐다.
기존 NexVirome 동등성, 자연 panel 성능, 전체 요구사항 완료를 의미하지 않는다.

## 1. 대조한 요구사항과 연구 목적

출처 식별자를 이후 작업에 유지한다.

- **S1**: 첫 첨부 자료(784c0dc8…): 보존 영역, taxonomic informativeness/VDS, masking의 손익,
  classification ambiguity와 assembly chimera의 관계, DB·read 조건 의존성, host/ERV, graph 연결 개선.
- **S2**: NMF 첨부 자료(2d1f8d9a…): 다중 샘플 coverage, composition·gene/domain·host 특징,
  단순 점수와 latent model 비교, 모호한 국소 graph에서 경로 선택, 계산량 제한.
- **S3**: 구조 분석 첨부 자료(8df99dfd…): sequence/fold/motif/context 기반 discovery,
  structural QC·binning, 후보 간 구조 비교, 실제 metagenome 전처리·복수 assembly·후속 재구성.
- **P1**: 사용자가 승인한 7단계 구현 계획: RefSeq snapshot, 자연 panel, ART,
  metaSPAdes/MEGAHIT, 경쟁 정렬 평가, 보수적 분할, 독립 평가 및 성공 기준.
- **U1**: 최근 명확화: 기존 NexVirome은 **reference에서 다른 taxonomy를 구별하는 데 도움이
  적은 영역을 임시 masking**한다. 공유서열 평가를 통해 그 좌표·범위·조건을 고도화하고,
  masking된 정보에 남아 있는 탐지·annotation·복원 가치를 확인한다.

연구는 세 축을 유지한다. (A) reference masking을 통한 분류 개선, (B) 잘못된 assembly 연결
재현·교정, (C) 보존 영역 재활용을 통한 바이러스 발견·annotation·복원이다. 같은 공유서열이
세 축에 주는 영향을 측정하되 분류 모호성과 assembly 오류가 같은 현상이라고 전제하지 않는다.
RdRP·RT·integrase 등의 명칭만으로 전체 유전자를 masking하지 않는다. 단백질 구조 보존도와
nucleotide 공유도 역시 별도 측정치다. 첨부 자료의 가중치·성능 수치·도구 순위는 검증된 기준으로 채택하지 않는다.

## 2. 요구사항별 구현 상태와 남은 산출물

| ID | 요구사항·출처 | 현재 확인된 상태 | 추가 산출물 |
|---|---|---|---|
| R01 | 버전 고정 RefSeq 원본 snapshot (P1) | 수집·대조 코드, 1개 accession smoke test | 전체 accession snapshot 실제 수집, annotation·taxonomy snapshot 확장, Linux 실행 |
| R02 | 자연 공유서열 panel (P1/S1) | 정확한 31-mer 후보 및 BLAST/flank 확인 | 위치별 공유 atlas, 분류 단계별 집계, circular·segment 그룹-aware panel |
| R03 | 기존 NexVirome mask 평가·환류 (U1/S1) | 미구현 | 원본/mask/정책 import·검증, 현행 baseline 재현, 새 정책 export·round-trip |
| R04 | 조건별 구별력/VDS (S1) | 미구현; 공유 k-mer 수는 대체 지표가 아님 | taxonomy rank·read 길이·library·DB별 구별력/오류율·표본수·불확실성 |
| R05 | masking의 이득과 손실 (U1/S1) | 미구현 | specificity/recall/미분류/정보 보존량 곡선, 영역별 유지·mask·재활용 근거 |
| R06 | 분류 모호성 대 assembly 오류 (S1) | assembly 평가만 있음 | 동일 원본 좌표로 결합한 2×2 분석, 조건별 관계 및 음성 결과 |
| R07 | 보존 영역의 생물학적 가치 (S1/S3/U1) | ORF·검색·구조 근거 일부 있음 | mask↔CDS/domain/motif 연결, 분류와 탐지·annotation 가치를 분리한 평가 |
| R08 | 숙주/ERV 및 비바이러스 대조군 (S1/S3) | 순수 바이러스 중심; 미구현 | host-only, ERV-only, exogenous-only, 혼합 및 진짜 integration 대조군 |
| R09 | 시퀀싱·strain·library 확장 (S1/P1) | PE150 및 coverage/ratio/seed 설계 | read 길이·insert·오류 profile·불균일 coverage·근연 strain·배경·library artifact 별도 실험 |
| R10 | 실제 다중 샘플 benchmark (S2) | coverage/NMF 실행 가능; seed 반복은 생물학적 변이 아님 | abundance를 독립적으로 바꾼 샘플 생성기, 상관·희소·비저차원 대조군 |
| R11 | 다중 특징 및 경로 점수 (S2/S3) | NMF는 coverage만, 단순 변화점 및 proposal gate | gene/domain·codon·host·구조 특징의 정규화/결측 처리, 단순 경로 점수와 ablation |
| R12 | taxonomic informativeness→graph (S1) | 미구현; 현재 교정은 reads/graph 중심 | 다중 좌표 대응 및 ambiguity annotation, reference 보조 arm과 reference-free arm 분리 |
| R13 | 모호한 reads·경로의 보존 (S1/S2) | 다중 정렬을 확정 지지에서 제외 | 가능한 원본/경로 집합·미해결 사유 보존, shared reads 활용과 확정 지지 분리 |
| R14 | 경로 선택·재구성 (S1/S2/S3) | 기존 contig 분할만; 제한된 merge–repeat–exit 탐지 | 보존된 경로의 순위·대안 출력, 검증 후 선택적 traversal/extension; 사라진 경로 복원은 후속 |
| R15 | 경쟁적인 키메라 평가 (P1) | consistent/chimeric/ambiguous/unaligned 및 gap-aware metrics | strain·segment·integration 대조, 식별 가능/불가능 복원 분리, circular 연속 복원 개선 |
| R16 | protein/domain과 motif 통합 (S1/S3) | MMseqs/Foldseek 및 Folddisco/FoldMason wrapper | CDS·profile/HMM annotation adapter, motif hit의 ORF·원본 좌표 및 curated label 연결 |
| R17 | discovery 네 축과 후보 family (S3) | 단순 후보 우선순위와 family Jaccard | sequence/global fold/local motif/context evidence 카드, motif-only·dark 구분, 구조 cluster 및 재현성 |
| R18 | 구조 계산 재사용·비용 제어 (S2/S3) | 검색 결과 SQLite cache, 국소 ORF 선택 | 독립 3Di/좌표 artifact cache, 후보 선별→고비용 분석, DB별 재검색에 예측 재사용 |
| R19 | geNomad/CheckV와 품질 통합 (S3/P1) | 실행 wrapper만 | native 결과 parser·ID mapping·통합 보고서; 도구 적용 범위와 unknown 보존 |
| R20 | 실제 데이터 전처리·assembly 운영 (S3) | reads→assembler wrapper | sample/host metadata, QC·선택적 host subtraction, per-sample 기본 흐름과 제거 reads 추적 |
| R21 | ensemble·dereplication·targeted coassembly (S3) | 두 assembler 독립 비교만 | 중복 제거 catalog·원본 대응, graph별 lineage 추적, 선택적 coassembly 비교 |
| R22 | 독립 평가와 누출 통제 (P1/S1–3) | accession 제외 및 기본 비교 | 근연 cluster/rank/segment 그룹 분리, annotation DB 누출 감사, 고정 정책·불확실성 평가 |
| R23 | 통합 보고 및 비교군 (P1/S2/S3/U1) | 기본 4군 metrics, 확장 5군은 주로 split count | masking·분류·교정·binning·discovery별 지표 및 각 교정군 matched random control |
| R24 | 자원·재현성·실패 재시작 (P1/S2) | manifest/Conda solve/65 tests/3 dry-run | 실제 Linux lock·실행, Python 포함 peak memory, 국소 mapping 재사용, stage 완료/실패 DAG 검증 |
| R25 | long-read/hybrid·별도 assembler (S2/S3/P1) | 의도적으로 후속 | 데이터와 기존 방법의 한계가 확인되면 비교; 독립 assembler 개발 의사결정 기록 |

현재 65개 테스트, 결과 파일 55개 동일성 및 workflow dry-run은 코드 검증이다.
자연 RefSeq의 오류 재현, 실제 masking 효용, 독립 panel의 개선 성능은 아직 입증되지 않았다.

## 3. P0 — 기존 NexVirome 계약과 평가 조건 고정 (R01/R03/R04/R22)

먼저 기존 NexVirome의 코드 revision, classifier 이름/버전/옵션, mask 생성 규칙,
원본 reference FASTA, 현행 masking 결과와 좌표, accession↔taxid 대응을 받아 재현한다.
현재 이 저장소에서 그 실제 구현·파일 형식은 확인되지 않았다. 추정한 방식으로 대체하지 않는다.

`ReferenceSnapshot`, `TaxonomySnapshot`, `MaskSet`, `ClassificationRun`을 독립 artifact로 정의한다.
taxonomy에는 accession.version·taxid·rank별 lineage, merged/deleted taxid 및 rank 누락 상태를 보존한다.
같은 taxon의 중복 reference 개수가 ambiguity를 부풀리지 않도록 accession 수와 distinct taxon 수를 구분한다.
표본 수와 DB 구성 편향도 함께 보고한다. 원본 FASTA는 불변이며 segment는 연결하지 않는다.

mask의 표준 내부 좌표는 원본 reference 기준 0-based half-open으로 한다. circular wrap은 두 구간과
동일 event ID로 표현한다. 원본 checksum·policy version·taxonomic rank·실험 조건·근거 ID를 필수로 둔다.
import 시 좌표 범위, accession version, 중첩 병합, 원본 서열을 검증하고 원래 형식으로 왕복 변환한다.

실제 classifier가 N masking/soft masking/seed exclusion을 어떻게 처리하는지 작은 fixture로 확인한다.
현재 `io.formats.fasta()`는 대문자로 읽고 기존 BLAST wrapper는 `soft_masking=false`이므로,
소문자 masking을 그대로 넣는 방식은 사용할 수 없다. mask-aware adapter를 별도로 구현한다.
길이를 바꾸는 삭제는 기본 정책으로 쓰지 않으며, 필요하면 명시적 좌표 대응표를 요구한다.

완료 기준: 현행 NexVirome 결과를 고정 fixture에서 재현하고, zero-mask/full-mask/overlap/circular/
잘못된 accession·checksum·길이/soft-mask 사례의 import·export·원본 복원이 검증된다.

## 4. P1 — 공유서열 atlas와 조건별 구별력 (R02/R04/R07/R09)

현재의 pair 후보 목록을 모든 관련 reference 위치로 확장한다. 정확한 canonical k-mer 후보 검색 뒤
경쟁 nucleotide 정렬을 사용하며, self hit·동일 taxon 공유·타 taxon 공유·반복·저복잡도를 구분한다.
정확한 31-mer가 없는 유사 구간도 놓칠 수 있으므로 별도 근사 후보 검색과 후보 recall 검증을 추가한다.
도구·k·정렬 기준을 기록하고 31-mer 하나를 보편적인 구별력 정의로 삼지 않는다.

RNA·DNA·retrovirus 및 taxonomy 거리별 strata를 별도로 선별·보고한다. 최초 상위 위험 pair 두 개가
모든 바이러스 유형을 대표한다고 보지 않는다. 비바이러스 경쟁 reference는 출처가 고정된 별도
decoy/background 집합으로 유지하며 원본 NCBI Virus RefSeq 집합과 섞어 동일 snapshot이라고 부르지 않는다.

구별력은 `D(region | rank, read_length, error_profile, PE_library, DB_snapshot, classifier)`처럼
조건을 명시한 측정으로 둔다. VDS라는 이름과 기존 수식이 있다면 P0에서 확인하고 호환 adapter를 둔다.
새 수식을 임의로 기존 VDS라고 부르지 않는다. 예측 ambiguity score와 실제 오분류율도 분리한다.

동일한 원본 reads/fragment를 경쟁 reference 전체에 평가해 correct-at-rank/wrong-at-rank/
ambiguous/unclassified를 보존한다. PE 단위와 개별 mate 결과를 구분하고 spanning unique flank가
공유 영역의 식별에 주는 이득을 계산한다. rank가 없거나 진짜 원본이 DB에 없는 사례를 별도로 집계한다.

산출물: `shared_regions.tsv`, `region_taxon_support.tsv`, `discriminability.tsv`, `domain_overlap.tsv`.
각 행에 original accession/좌표, 경쟁 accession·taxon, identity·길이, 반복 유형, 평가 조건,
관측 수와 신뢰구간을 포함한다. 영역의 문제성이 단백질 이름에서 추론된 것인지 실제 reads에서 측정된
것인지 구분한다. full ORF/domain을 구조 분석에 사용하고 짧은 nucleotide window별 구조예측은 피한다.

완료 기준: 동일 species/서로 다른 species/같은 genus/다른 family, reverse complement,
반복·저복잡도, 원형·segment, shared-only와 unique-flank, read 길이 변화 및 DB 추가 시 기대 변화 검증.

## 5. P2 — masking 최적화·재활용·기존 NexVirome 환류 (R03/R05/R07)

같은 원본 read set과 classifier 설정으로 다음 군을 비교한다.

| 군 | 목적 |
|---|---|
| 원본 reference | masking 없는 분류·탐지 기준 |
| 기존 NexVirome mask | 반드시 재현할 실제 baseline |
| nucleotide 공유도 기반 정책 | 구조 정보 없이 가능한 개선 |
| rank/library/DB 조건별 정책 | 실제 구별력에 맞춘 mask 범위·강도 |
| 위 정책 + 선택적 보존 영역 재활용 | 고유 근거로 후보를 좁힌 뒤 탐지·annotation·복원에 활용 |
| reference별 같은 양/구간 길이 분포의 random mask | 단순 정보 제거량 효과를 통제 |

masking은 reference의 분류용 view에 적용한다. **모든 군의 simulation은 동일한 unmasked 원본에서
만들고, raw reads를 masking 정책 때문에 삭제하거나 바꾸지 않는다.** masked reference에서 reads를
생성해 sensitivity 손실을 감추지 않는다. assembly 연결 개선 실험에서도 원래 reads를 유지한다.

정책은 우선 설명 가능한 threshold로 구현한다. mask 양·남은 informative breadth·recall 제약하에서
오분류 감소를 비교하고, 개발 panel에서 선택한 임계값과 허용 손실을 평가 전에 고정한다.
0 관측 영역은 근거 부족으로 남기며 과도한 masking으로 FP만 줄인 정책을 성공으로 보지 않는다.
유전자/domain 경계와 겹침을 보고하되 경계 전체를 자동 확장 masking하지 않는다.

재활용에서는 고유 영역이 지지하는 taxon 후보와 보존 영역의 호환성을 확인한다. 보존 영역만으로
하위 taxon을 확정하지 않는다. 후보 밖의 경쟁 설명도 audit에 유지해 후보 축소가 강제 오분류를
만들지 않는지 검사한다. 분류 가중치·탐지 근거·annotation 근거는 각각 기록한다.

지표: rank별 precision/recall, FP율(분모 명시), 미분류·모호 비율, macro/micro 성능,
masked fraction, 남은 식별 가능 breadth, 추가 탐지·ORF/domain 회수량, 원본별 recovery.
전체 입력을 분모에 포함하며 할당된 reads만으로 sensitivity를 계산하지 않는다.
일반화는 reference/pair/sample 단위로 평가하고 reads 수를 독립 생물학적 반복 수로 취급하지 않는다.

산출물: `mask_proposals.tsv`, `policy.json`, `mask_diff.tsv`, `mask_benefit.tsv`, NexVirome 호환 export.
기존 mask는 덮어쓰지 않고 versioned candidate를 생성한다. export→import round-trip과 기존 classifier
실행으로 검증한다. 실험 DB·classifier 버전이 다르면 기존 score 재사용을 차단하거나 명시적으로 재평가한다.

완료 기준: 현행 대비 이득/손실 곡선, 독립 panel 결과, 정보 보존 조건, 실패·음성 결과와 정책 변경 이유가
포함된 review 가능한 정책 패키지를 만든다. masking의 수치 성공 목표는 현행 baseline 측정 뒤 고정한다.

## 6. P3 — 분류 오류와 assembly 오류의 연결 검증 (R06/R08/R09/R12/R15)

동일 원본·동일 mixture에서 분류 평가와 두 assembler의 평가를 실행한다. contig↔원본 정렬의
다중 대응과 breakpoint interval을 보존해 공유 atlas에 투영한다. 고유하지 않은 위치는 임의의
원본 좌표로 확정하지 않는다. mask 구간·분류 오류·graph 분기·확정 잘못된 연결을 하나의 표로 연결한다.

분류 모호성 유/무 × assembly 오류 유/무를 조건별로 보고한다. 관측된 overlap·enrichment와
불확실성을 측정하고 길이·반복성·coverage·유사도 차이를 통제한다. 관련성이 낮거나 한 축에서만
이득이 있어도 그대로 결과로 남긴다. 모든 mask 구간을 chimera-prone 구간으로 간주하지 않는다.

host-only/ERV-only/exogenous-virus-only/혼합 실험을 별도 strata로 추가한다. human·mouse·bat는
원자료에 나온 일반화 후보이며 실제 숙주 metadata와 버전 고정 reference를 확인한 뒤 선택한다.
human ERV는 HERV, 다른 숙주는 해당 ERV로 표기한다. 실제 provirus/integration·재조합 서열을
단일 원본 정답으로 넣어 생물학적 host-virus 경계를 assembly 오류로 오인하지 않는지 확인한다.
숙주 제거 전후 비교에서는 제거된 reads와 false viral loss를 추적한다.

P1 기본 PE150/300±30, ART profile, 50×·1:1·3 seeds 및 10/50/100×·1:1/1:10/1:100 확장은 유지한다.
추가 read 길이(예: 100/150/250), insert 조건, 오류 profile, strain 혼합, 불균일 coverage, index hopping
등은 지원 여부를 확인한 명시적 실험으로 분리한다. ART 오류와 별도 오류율을 중복 적용하지 않는다.
assembler의 실제 k 정책은 구조화해 수집하고 동일 숫자 k를 동일 알고리즘 조건이라고 해석하지 않는다.

독립성은 accession뿐 아니라 근연 sequence cluster·genome/segment 묶음으로 관리한다.
DB-version 비교는 동일 read set으로 수행하며, closed-set·held-out sequence·held-out taxon 실험을
구분한다. source truth는 simulation/evaluation에만 둔다. reference 보조 교정군은 허용된 공개 DB만
받고, reference-free baseline과 분리해 보고한다. structure/sequence annotation DB의 held-out 중복도 감사한다.

## 7. P4 — 다중 샘플·graph 연결 판정 고도화 (R10–R14/R18)

서로 다른 생물학적 abundance를 갖는 sample design을 추가한다. seed 변경은 technical replicate다.
같이 움직이는 원본, 독립 원본, sparse/zero coverage, 같은 비율로 움직여 분리가 불가능한 원본,
low-rank 가정이 약한 경우를 모두 포함한다. 공통 경쟁 target과 library normalization을 고정한다.

현재 coverage-only NMF와 cosine 기반 변화점, read-pair baseline을 보존한다. 추가 특징은 coverage,
k-mer/GC/codon usage, gene/domain, host homology, 구조 profile로 분리하고 block별 scaling·결측·음수
처리를 명시한다. 현재의 빈 구조 근거는 반대 근거가 아니다. 단순 weighted path score와 각 특징
ablation을 먼저 비교하고 NMF 결합 모델은 재구성 오차만이 아니라 실제 held-out 오류/복원 성능으로 선택한다.

현재 proposal gate는 이미 read 근거가 있는 split을 제한할 뿐, 대안 경로 자체를 선택하지 않는다.
다음 단계에서는 upstream/shared/downstream·대체 경로별 feature와 rank별 ambiguity를 함께 출력한다.
불확실한 reads의 가능한 경로 집합을 보존하고 확정 fragment 지지와 구분한다. 근거가 부족하면
`unresolved`를 유지한다. 같은 원본에 호환되는 양쪽 flank는 공유 영역을 유지하는 근거로 사용한다.

안전한 split을 기준선으로 유지하면서, graph에 실제로 남아 있는 경로만 선택·재구성하는 후속 모드를
구현한다. 새로운 염기를 보충하지 않고 경로 방향·overlap·좌표 대응 및 read 검증을 요구한다.
복잡한 bubble/cycle/긴 반복/누락 경로는 지원 범위와 실패 사유를 명시한다. 단순화로 이미 사라진
경로 복원이나 독립 assembler는 별도 단계이며 지금의 분할 prototype 완료와 혼동하지 않는다.

계산은 모호한 neighborhood에 집중한다. 현재 후보마다 full read library를 재검색하는 비용을
줄이기 위해 공통 mapping/index와 read 추출을 재사용하되, 경쟁 경로 누락으로 specificity가 과장되지
않도록 full-search 대조를 둔다. 행렬/후보 수 제한, 배치·디스크 인덱스 및 비용 측정을 추가한다.

## 8. P5 — masking 영역 재활용과 구조 discovery 완성 (R07/R16–R21)

현재 six-frame ORF는 후보 생성기다. RefSeq CDS와 실제 viral gene caller의 검증된 결과를
좌표·strand·genetic code·partial 상태와 함께 import한다. profile/HMM 결과 adapter를 추가하여
full protein/domain 단위의 sequence·global fold·local motif·genome context를 구분한다.
RdRP와 retroviral RT/Pol 등을 같은 marker로 취급하지 않는다. motif의 원자 좌표와 잔기 번호는
protein/ORF/reference nucleotide 위치로 명시적으로 대응한다.

Folddisco는 현재 raw motif 결과, FoldMason은 native alignment/report wrapper다. 결과 parser,
confidence·negative control·curated label provenance를 구현한 뒤 근거에 통합한다. 검색을 안 한 것,
실패한 것, 신뢰도가 낮은 예측, 진짜 no-hit을 구분한다. predicted 3Di를 원자 좌표처럼 쓰지 않는다.

후보 유형은 known sequence, remote fold relative, motif-only candidate, uncharacterized/dark로
나누고 gene order/synteny·host similarity·hallmark 조합·샘플 간 재현성을 함께 보고한다.
순위 점수는 확률이나 신규 species의 증명이 아니다. 유사 구조를 가진 후보들을 비교·cluster하고
진화적 해석은 독립 alignment·context 근거가 있을 때만 한다. no-hit만으로 바이러스라고 선언하지 않는다.

ORF별 3Di/좌표 예측 artifact를 서열·모델·버전·설정·품질 지표로 캐시하여 DB가 바뀌어도 재사용한다.
가벼운 sequence/profile 단계→의심 후보 선정→3Di→필요 후보의 실제 구조 및 motif 순으로 비용을 관리한다.
현재 검색 결과 cache와 예측 cache를 구분한다. Phold는 별도 backend 후보이고 현재 경로는 Foldseek/ProstT5다.

geNomad/CheckV native 결과를 공통 ID/좌표와 연결하고, 바이러스 유형에 맞지 않는 completeness 값은
unknown으로 남긴다. 구조 family는 taxonomy나 genome membership과 동일시하지 않는다.
binning은 fragment association이며 정렬 순서나 genome 연결을 의미하지 않는다. simulation truth로
bin purity/completeness·잘못 섞인 source를 평가하고 구조 가중치 없는 대조군을 유지한다.

실제 데이터에는 sample/host manifest, adapter/quality QC, 선택적 host subtraction, per-sample assembly,
read 재정렬 흐름을 추가한다. 제거 reads·contigs와 이유를 남기고 모든 sample을 human으로 가정하지 않는다.
세균·진균·고세균·unknown 서열을 숙주 제거라는 이유만으로 일괄 제거하지 않는다. 특히 prophage나
숙주와 공유하는 viral gene을 잃는 경우를 negative/retention control로 검사하고 제거 정책을 명시한다.
여러 assembler의 결과를 dereplication한 catalog는 원본 contig·assembler·graph 경로를 역추적할 수 있어야 한다.
두 assembler의 일치를 정답으로 간주하지 않으며, targeted coassembly는 strain 혼합 위험과 함께 비교한다.

## 9. P6 — 독립 비교·통합 보고·실제 실행 (R15/R19/R22–R24)

기존 MEGAHIT/metaSPAdes/corrected/matched random split 비교를 유지한다. 확장 baseline/simple/NMF/
structure/combined 모두에 실제 evaluate와 자기 split 수에 맞춘 random control을 연결한다.
현재 확장 workflow의 split count만으로 개선을 주장하지 않는다. masking 비교군은 별도 표로 유지한다.

보고서는 (1) masking·분류, (2) 키메라·연결 교정, (3) 식별 가능한/호환되는 원본 recovery,
(4) annotation/discovery/binning, (5) runtime·memory·저장공간·cache 효과를 분리한다.
무필터/500 bp 결과와 모든 짧은 조각을 보존한다. 원형 접합부를 통과하는 올바른 연속 복원과
segment별/그룹별 복원을 구분하고, scaffolds는 contig와 별도 평가한다.

개발 panel에서 정책·threshold를 고정하고 독립 panel에서 confidence interval과 조건별 실패를 보고한다.
assembly 목표는 확정 잘못된 연결 50% 이상 감소, recovery 감소 2 percentage points 이내를 유지한다.
baseline 오류가 0이면 감소율을 정의하지 않는다. 자연 pilot→정해진 확장까지 오류가 없으면
관찰 범위에서 재현되지 않았다고 보고한다. masking 연구는 계속할 수 있으나 assembly 개선 주장으로 바꾸지 않는다.

Linux에서 실제 Conda 설치·doctor·최소 end-to-end 실행 뒤 정확한 패키지 lock, pip 설치와 code revision,
DB·모델 checksum을 기록한다. Python 부모와 subprocess를 포함한 peak memory, 실제 k·도구 버전·seed·실패
로그를 수집한다. 성공 산출물만 DAG 완료로 인정하는 failure/restart fixture를 넣고, 실패한 stage는
기존 결과를 조용히 덮어쓰지 않는 재시도 절차를 유지한다.

## 10. 모듈·CLI·workflow 추가 계획

아래는 **새로 구현할 인터페이스 설계**이며 현재 실행 가능한 CLI로 표시하지 않는다.

| 추가 모듈 | 책임 | 계획한 CLI/산출물 |
|---|---|---|
| `taxonomy/` | rank lineage snapshot 및 taxid 정규화 | `taxonomy build/validate` |
| `sharing/` | 위치별 공유 atlas와 조건별 구별력 | `sharing build/evaluate` |
| `masking/` | import, 정책 후보, masking view, 버전 비교 및 export | `mask import/propose/evaluate/export` |
| `classification/` | 실제 NexVirome classifier adapter, rank별 경쟁 근거·성능 | `classify`, `classification evaluate` |
| `benchmark/` 확장 | mask/host/multisample/DB drift/held-out 설계 | 생성된 실험 manifest 및 정답 별도 저장 |
| `evidence/` 확장 | mask↔ORF/domain/motif↔graph 좌표 및 근거 카드 | `evidence link` |
| `correction/` 확장 | reference 보조 경로 점수와 보존 경로 재구성 | 별도 opt-in 경로 모드 및 좌표 대응 |
| `discovery/`·`catalog/` | 후보 유형, structure cluster, dereplication·bin 검증 | catalog·근거 카드·오류/복원 보고 |
| `application/`·workflows | 동일 use case로 모든 비교군 실행 | `masking.smk`, `joint_benchmark.smk`, 실데이터 workflow |

`io/domain/runtime`는 기존 리팩터링 구조를 재사용한다. `MaskPolicy`, `ClassifierBackend`,
`TaxonomyResolver`는 교체 가능한 계약으로 두며 파일 parsing과 점수 계산을 분리한다.
새 도구가 필요한 경우 해당 단계 구현 시 공식 CLI·버전·Conda 설치를 검증한다. 자료에 등장한
DIAMOND/HMMER/Phold/COBRA 등의 모든 대안을 최초 필수 의존성으로 동시에 추가하지 않는다.

## 11. 후속·조건부 항목과 착수 순서

Long-read/hybrid(HiFi/ONT와 관련 assembler), COBRA류 extension 비교, MEGAHIT graph–contig
대응, viral-specific SPAdes 모드, RNA-seq 특성, sequence embedding/VAE/GNN/transformer,
strain/haplotype 동시 복원, unseen-virus 예측, 별도 assembler는 누락하지 않고 후속 목록에 둔다.
실제 데이터 제공, 단순 모델의 성능, 남아 있는 오류 유형과 자원 측정을 근거로 착수한다.
도구의 최신 순위나 자료 속 성능 배수는 계획의 전제가 아니다.

착수 순서: **P0 계약·baseline → P1 공유 atlas → P2 mask 평가·export → P3 분류/assembly
공동 검증 → P4 경로·다중 샘플 개선 → P5 discovery·실데이터 통합 → P6 독립 보고**.
P6의 환경 검증과 평가 fixture는 각 단계와 함께 준비한다. 각 작업은 R ID, 입력/출력 계약,
완료 테스트, 실행 환경, 실제 실험 상태를 기록해 다음 대화에서 다시 범위가 빠지지 않게 한다.

구현 착수에 필요한 미확인 입력: 기존 NexVirome 코드·mask 형식·classifier 설정·reference 버전,
실제 sample/host/library metadata와 reads 위치, Linux 실행 환경, curated domain/structure labels.
이 입력의 부재는 계획에 기록하며 임의의 현행 masking 구현이나 숙주를 만들어 넣지 않는다.
