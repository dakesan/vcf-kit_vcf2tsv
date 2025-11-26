# vcf2tsv

Convert VCF files to TSV format.

## Requirements

- Python >= 3.8
- bcftools (must be installed and available in PATH)

## Installation

```bash
uv pip install .
```

Or with uvx:

```bash
uvx --from . vcf2tsv
```

## Usage

```bash
# Wide format (default)
vcf2tsv input.vcf

# Long format
vcf2tsv input.vcf --format long

# With header
vcf2tsv input.vcf --header

# Expand snpEff ANN annotations
vcf2tsv input.vcf --header --ann

# Show help
vcf2tsv --help
```

## Options

- `--format, -f`: Output format (`wide` or `long`). Default: `wide`
- `--header, -h`: Print header line
- `--ann, -a`: Expand snpEff ANN annotations

## License

MIT
