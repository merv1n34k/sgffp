.PHONY: setup dev build test test-all lint fmt clean

setup:
	uv sync --all-extras

dev:
	uv run sff --help

build:
	uv build

test:
	uv run pytest tests/ -v

test-all:
	uv run pytest tests/ -v --cov=sgffp

lint:
	uv run ruff check .

fmt:
	uv run ruff format .

clean:
	rm -rf dist/ .pytest_cache/ __pycache__/ .ruff_cache/ .coverage htmlcov/ src/*.egg-info
