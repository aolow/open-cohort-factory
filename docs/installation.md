# Installation

## Prerequisites

- macOS, Linux, or Windows
- Git
- Python 3.11 or newer, either already installed or managed automatically by uv
- Internet access when querying public data sources

No database, cloud account, or GDC authentication token is required for the current public-metadata
workflow.

## Install uv

On macOS with Homebrew:

```bash
brew install uv
```

On macOS or Linux without Homebrew:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

On Windows PowerShell:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

## Clone and install

```bash
git clone https://github.com/aolow/open-cohort-factory.git
cd open-cohort-factory
uv sync --all-extras --dev
```

The repository is currently private. Anyone cloning it must first be granted access by its owner
and authenticate with GitHub.

## Verify the installation

```bash
uv run cohort-factory --help
uv run cohort-factory validate examples/luad_older_adults.yaml
uv run pytest
```

## Build the documentation locally

```bash
uv run mkdocs serve
```

Open the local address printed by MkDocs. To perform the same strict documentation check used in
development:

```bash
uv run mkdocs build --strict
```

