from pathlib import Path


def test_project_layout() -> None:
    root = Path(__file__).parents[1]
    for relative in ("configs/reference.yaml", "configs/simulation.yaml", "docs/study_design.md"):
        assert (root / relative).exists()

