# Installation

Follow the [basic requirements](#basic-requirements) below to install and run outward assembly.

By default, the entire pipeline runs on your local machine. For very large datasets, there's an option to run the read search step on AWS Batch, trading startup time for potentially very wide parallelization. See the [usage docs](./usage.md#optional-aws-batch-for-read-search) to learn when this might be helpful, and [AWS Batch setup](#optional-aws-batch-for-read-search) for additional installation steps.

## Basic requirements

Outward assembly only supports Linux.

### Python dependencies

Python dependencies are managed with [uv](https://docs.astral.sh/uv/). [Install](https://docs.astral.sh/uv/getting-started/installation/) uv if you don't already have it. Then from the base repository directory:

```bash
# Install and sync Python dependencies with uv
uv sync --extra dev

# Run Python commands with uv
uv run your_script.py
uv run pytest
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

### AWS CLI and credentials

Outward assembly streams reads from S3, which requires the [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html) and credentials with permission to list and read your S3 objects.

Configure credentials by creating a file at `~/.aws/credentials`:
```
[default]
aws_access_key_id = <ACCESS_KEY_ID>
aws_secret_access_key = <SECRET_ACCESS_KEY>
```

See the [AWS CLI documentation](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html) for more configuration options.

## (Optional) AWS Batch for read search

Outward assembly uses Nextflow as a (somewhat heavy) frontend to AWS Batch, so running read searches on Batch requires installing Nextflow and setting up AWS Batch infrastructure.

### Install Nextflow and Docker

Install the following:
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Docker](https://docs.docker.com/engine/install/)

Make sure your user is configured to use Docker by adding yourself to the `docker` group:

```bash
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

> [!TIP]
> If you encounter `AccessDenied` errors when running Nextflow, you may need to export your AWS credentials as environment variables:
>
> ```bash
> eval "$(aws configure export-credentials --format env)"
> ```

### Seqera Platform access token

Running read search on AWS Batch requires a Seqera Platform [access token](https://docs.seqera.io/platform-cloud/api/overview). See your tokens or create a new one in the Seqera [console](https://cloud.seqera.io/tokens).

### Set up AWS Batch

You'll need a Batch job queue and compute environment. You can use existing Batch infrastructure or follow the steps below to set it up via the AWS console.

**If using existing infrastructure:** Outward assembly creates many short-lived Batch jobs (hundreds of thousands of jobs, each running for a few seconds). This workload doesn't need large or high-performance EBS volumes—free-tier gp3 settings (3000 IOPS, 125 MB/s) are sufficient. Avoid over-provisioning storage to control costs.

#### 1. Check your AWS permissions

The most common problems when setting up Batch arise from insufficient AWS permissions. Make sure you have:

- **AWS Batch:** Permissions to create and manage compute environments and job queues (e.g. `AWSBatchFullAccess` policy).
- **EC2:** Permissions to create launch templates and compute environments (e.g. `AmazonEC2FullAccess` policy).
- **Instance role with S3 access:** Batch containers inherit permissions from the compute environment's instance role. This role needs S3 read access to your data buckets plus ECS permissions (e.g. `AmazonS3ReadOnlyAccess` + `AmazonEC2ContainerServiceforEC2Role`, or ask your administrator for an appropriate role).

#### 2. Create an EC2 launch template

Navigate to **EC2 > Launch templates** in the AWS console:

1. Click "Create launch template".
2. Enter a name (e.g. `outward-assembly-batch-template`).
3. Check the box under "Auto scaling guidance".
4. Under "Application and OS Images", click "Browse more AMIs", search for "Amazon ECS-Optimized Amazon Linux 2023 (AL2023) x86_64 AMI", and select it from AWS Marketplace AMIs. Subscribe if prompted (it's free).
5. Select an instance type (e.g. `c6a.xlarge`). This is just a default—Batch will provision different types based on the allowed instance types in your compute environment.
6. Under "Key pair (login)", select "Don't include in launch template".
7. Under "Network settings", select "Create security group" with default settings.
8. Under "Storage (volumes)", expand the default volume and configure:
   - **Size:** 30 GiB (enough for the AMI plus working space)
   - **Volume type:** gp3
   - **IOPS:** 3000 (free tier)
   - **Throughput:** 125 MB/s (free tier)
9. Add tags for cost tracking (e.g. an "Instances" tag and a "Volumes" tag).
10. Click "Create launch template".

#### 3. Create a Batch compute environment

Navigate to **Batch > Compute environments** in the AWS console:

1. Click "Create".
2. Select "Amazon Elastic Compute Cloud (Amazon EC2)" and "Managed" orchestration.
3. Enter an environment name (e.g. `outward-assembly-compute`).
4. Set up roles:
   - **Service role:** `AWSServiceRoleForBatch`
   - **Instance role:** The role from step 1 with appropriate permissions
5. Under "EC2 tags", add a tag to identify your usage (e.g. `outward-assembly`).
6. Click "Next". Under "Instance configuration":
   - Enable "Use EC2 Spot instances" (recommended for cost savings).
   - Set "Maximum % on-demand price" to 100.
   - Set vCPUs: Minimum 0, Desired 0, Maximum 1024.
   - Under "Allowed instance types", select compute-optimized families like `c6a`, `c6i`, `c7a` (read search benefits from fast CPUs and network, not large memory).
   - Under "Additional configuration", select your launch template.
7. Click "Next" twice (accept network defaults), then "Create compute environment".

#### 4. Create a Batch job queue

Navigate to **Batch > Job queues** in the AWS console:

1. Click "Create".
2. Select "Amazon Elastic Compute Cloud (Amazon EC2)".
3. Enter a queue name (e.g. `outward-assembly-queue`).
4. Select your compute environment from the dropdown.
5. Click "Create job queue".
