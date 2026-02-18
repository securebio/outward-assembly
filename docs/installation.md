# Installation

To run this software, you need to follow the instructions in [basic requirements](#basic-requirements). This will allow you to run the outward assembly software in the *local* profile (we expect the *local* profile will work for most use-cases).

Briefly, the outward assembly software has two profiles (read our [usage docs](./usage.md#profiles) to learn more about when to use either profile):
* *Local*: Run outward assembly using your local machine.
* *Batch*: Run outward assembly using AWS Batch (via Nextflow) for distributed read searching; other algorithm steps run locally.

If you want to use the *batch* profile, you'll also need to follow the [instructions below to install the required software](#optional-batch-profile).

## Basic requirements

### Python dependencies

Python dependencies are managed with [uv](https://docs.astral.sh/uv/). From the base repository directory:

```bash
# Install and sync Python dependencies with uv
uv sync --extra dev

# Run Python commands with uv
uv run your_script.py
uv run pytest

# Or activate the virtual environment
source .venv/bin/activate  # Linux/macOS
# .venv\Scripts\activate     # Windows
```

### Bioinformatics tools

Bioinformatics tools are managed via conda using the **tools-only** environment (`oa_tools_env.yml`) that contains just the bioinformatics tools (MEGAHIT, fastp, KMC) without Python dependencies.

**Additionally, you need to install Nucleaze** (a Rust-based k-mer filtering tool). See [nucleaze_migration.md](./nucleaze_migration.md) for installation instructions, or briefly:

```bash
# Requires Rust 1.88+
git clone https://github.com/jackdougle/nucleaze.git
cd nucleaze
cargo install --path .
```

Install the conda tools environment:

```bash
mamba env create -n oa-tools -f oa_tools_env.yml --channel-priority flexible
```

**Before running the pipeline, activate the tools environment:**

```bash
mamba activate oa-tools
```

Then run Python commands via uv (while the tools environment is activated):
```bash
uv run your_script.py
```

**Note:** You need both uv (for Python packages) and mamba/conda (for bioinformatics tools like MEGAHIT, KMC, etc.), plus Nucleaze installed via cargo.

## (Optional) Batch profile

Using the batch profile requires doing two steps:

1. Install Nextflow and the AWS CLI
2. Configuring AWS Batch

### Nextflow installation

#### 1. Install Nextflow 
To install Nextflow, you can follow the instructions [here](https://www.nextflow.io/docs/latest/getstarted.html)

#### 2. Install AWS CLI
To install the AWS CLI, you can follow the instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

#### 3. Install Docker
To install Docker, you can follow the instructions [here](https://docs.docker.com/engine/install/).

#### 4. Configure AWS & Docker

[Configure AWS access](https://www.nextflow.io/docs/latest/aws.html) by creating a file at `~/.aws/config` **or** `~/.aws/credentials`, specifying your access key ID and secret access key, e.g.

`~/.aws/config`:
```
[default]
region = us-east-1
output = table
tcp_keepalive = true
aws_access_key_id = <ACCESS_KEY_ID>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

`~/.aws/credentials`:
```
[default]
aws_access_key_id = <ACCESS_KEY_ID>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

> [!TIP]
> If you encounter `AccessDenied` errors after doing this, you may also need to export these keys as environment variables before running Nextflow:
>
> ```
> eval "$(aws configure export-credentials --format env)"
> ```

Next, you need to make sure your user is configured to use Docker. To do this, create the `docker` user group and add your current user to it:

```
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

### Setup AWS Batch

To setup AWS Batch, we recommend that you follow the latest [^1] instructions in the [MGS Workflow](https://github.com/naobservatory/mgs-workflow/blob/master/docs/batch.md) [^1], _except_ substitute the following storage configuration in your EC2 launch template:
* Volume type: gp3
* Capacity: 5GB plus whatever space is required for the snapshot used by your launch template
* Throughput: 125MB/s (free tier)
* IOPS: 3000 (free tier)

Though each Batch job runs for just a few seconds, running outward assembly at scale may create hundreds of thousands of Batch jobs. Thus, to control costs, it's important not to provision unnecessarily large EBS volumes. 

[^1]: In the situation that the instructions change, [this link](https://github.com/naobservatory/mgs-workflow/blob/9fe05a5ca9ce7cbc886927788f22c71ff9f26443/docs/batch.md) should be used as reference for the version of the docs used at the time.

### Seqera Tower Credentials
Using the Batch profile requires a Seqera Tower [access token](https://docs.seqera.io/platform-cloud/api/overview). See your tokens or create a new one in the Seqera [console](https://cloud.seqera.io/tokens).

 
