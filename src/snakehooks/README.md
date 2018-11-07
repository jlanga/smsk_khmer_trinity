[![Build Status](https://travis-ci.org/jlanga/snakehooks.svg?branch=master)](https://travis-ci.org/jlanga/snakehooks)

# snakehooks

git hooks for snakemake pipelines

## Getting started

Writing maintainable pipelines requires automatic checks to ensure:

- That the syntax is OK.
- That the scripts and source code adheres to some type of language conventions
    to prevent smelly code.
- That the documente

### Prerequisites

- `snakemake`,
- `shellcheck`,
- `yamllint>=1.11.1`: the `--strict` option happens here
- `pylint`

### Installing

Link or move `pre-commit` to `.git/hooks/`, or run
`./src/snakehooks/link_hooks` supposing that this repo is installed in `src/`.

### Running it

```sh
git commit
```

That's all. The tests will be run in the following order:

- `snakemake clean` to remove previous builds,
- `snakemake --dryrun` to check the snakefiles,
- `pylint` over the Snakefiles and the subworkflows in `src/snakefiles/`,
- `pylint` over `*.py` scripts that are not in `src/snakefiles/`,
- `yamllint` over `*.{yaml, yml}` files, and
- `snakemake` to run the pipeline from beginning to end.

### Contributing

Do you want to propose another tool? To modify an existing one? Fork and make
a PR!
