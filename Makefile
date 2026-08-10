.PHONY: install test lint format typecheck check

install:
	pip install -e ".[dev]"

test:
	pytest

lint:
	ruff check .

format:
	black .

typecheck:
	mypy src

check:
	ruff check .
	black --check .
	mypy src
	pytest
