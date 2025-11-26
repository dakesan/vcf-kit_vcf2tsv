# vcf2tsv

Convert VCF files to TSV format.

This is a minimal standalone package extracted from [VCF-kit](https://github.com/AndersenLab/VCF-kit).

## Requirements

- Python >= 3.8
- bcftools (must be installed and available in PATH)

## Installation

```bash
pip install git+https://github.com/dakesan/vcf-kit_vcf2tsv.git
```

Or clone and install locally:

```bash
git clone https://github.com/dakesan/vcf-kit_vcf2tsv.git
cd vcf-kit_vcf2tsv
pip install .
```

## Usage

```bash
# Wide format (default)
vcf2tsv input.vcf

# With header
vcf2tsv input.vcf --header

# Long format (one sample per row)
vcf2tsv input.vcf --format long --header

# Expand snpEff ANN annotations
vcf2tsv input.vcf --header --ann

# Output to file
vcf2tsv input.vcf --header > output.tsv

# Show help
vcf2tsv --help
```

## Options

| Option | Short | Description |
|--------|-------|-------------|
| `--format` | `-f` | Output format (`wide` or `long`). Default: `wide` |
| `--header` | `-h` | Print header line |
| `--ann` | `-a` | Expand snpEff ANN annotations |

## Output Formats

- `wide`: One variant per row. Samples are expanded as columns.
- `long`: One variant x sample per row. Samples are expanded as rows.

## License

MIT
