FROM julia:1.11.4

WORKDIR /app

# Copy environment definition first (so installs are cached)
COPY Project.toml .
# COPY Manifest.toml .  # optional: include if you’ve pre-generated it locally

# Install & precompile environment
RUN julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

# Copy all scripts and utils
COPY *.jl .

# Default: open Julia REPL in this environment
CMD ["julia"]
