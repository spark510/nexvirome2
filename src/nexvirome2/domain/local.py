from dataclasses import dataclass

@dataclass
class LocalPaths:
    templates: dict[str, str]
    bounds: dict[str, tuple[int, int]]
    paths: list[list[str]]
