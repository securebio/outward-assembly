# Developer notes

This is a collection of notes for developers working on the outward assembly pipeline.

## Nomenclature

The following nomenclature is used in the outward assembly pipeline:
* outer iteration: An outer iteration is a single run of the outward assembly pipeline and will only be referenced when we are running outward assembly in *automated* mode.
* inner iteration: An inner iteration refers to the amount of iterations it takes for one run of the outward assembly pipeline to complete.

Hence, *normal* mode will only have one outer iteration and a variable number of inner iterations, while *automated* mode will potentially have multiple outer iterations and different amounts of inner iterations for each outer iteration.

## Code notes
* The naive algorithm used to construct the contig overlap graph requires checking all pairs of contigs for overlap and thus has runtime `O(n^2)` in `n` the number of contigs. Typically this step is super fast because we have few contigs, but if you're finding outward assembly stalling with a single Python process utilizing ~100% of a CPU core, this is a plausible candidate.
    * One way to accidentally get lots of contigs: use a too-small `read_subset_k` argument when calling `outward_assembly`, which can result in BBDuk selecting lots of reads unrelated to the seed sequence.
* For code formatting, we use [Black](https://github.com/psf/black) and [isort](https://github.com/PyCQA/isort).

## Testing

We use the [pytest](https://docs.pytest.org/en/stable/) testing framework. Tests should be annotated with the markers defined in [pyproject.toml](../pyproject.toml).

You can run the full test suite with `pytest` in the base repo dir, or `pytest -s` to see the real time logs. The full suite includes integration and end-to-end tests, and can be quite slow. You may find it useful to run just the tests with certain markers, e.g. `pytest -m unit` to run just the unit tests.

Some integration and end-to-end tests run in both local and Batch mode, hence you'll have to have gone through the [AWS Batch setup](../docs/installation.md#optional-aws-batch-for-read-search) to run the full test suite. 

To run the tests that use Batch and AWS, copy `.env.example`  to `.env` and fill in your compute resources and credentials. (Alternatively, you can directly export the variables defined in `.env` as environment variables.)

## Improvements
* [Google doc](https://docs.google.com/document/d/1AiQUWMNUhbwYZBqLleZ1K-XXnND84z0tURidb2OD8sw/edit?tab=t.0) with some ideas for performance and sensitivity enhancements.
* Richer test coverage; existing test coverage is pretty modest.
* Set up GitHub actions testing
* Rearranging code for greater clarity and shorter files
* Better logging:
    * More comprehensive
    * Principled usage of log levels
    * Optionally not discarding all logs from parallel processes (currently written to `/dev/null`)
* Refactor use of working directory so that every iteration of outward assembly has its own subdirectory and we're not updating files like `current_contigs.fasta` in a stateful way.
