# CLAUDE.md - Guidelines for mc-slim repository

## Build & Run Commands
- Build Docker container: `docker build -t my_openmc_custom .`
- Run container: `docker run -it --name=my_openmc_custom -v src:/root/src my_openmc_custom`
- Start container: `docker start -i my_openmc_custom`
- Run Python script: `python3 filename.py` (inside container)

## Code Style Guidelines
- **Imports**: Standard library first, then third-party, then local modules
- **Formatting**: 4-space indentation, 100-character line limit
- **Types**: Use type hints where appropriate for function parameters and returns
- **Naming**: snake_case for functions/variables, CamelCase for classes
- **Comments**: Section headers with `###` comment blocks
- **Error handling**: Use try/finally blocks for resource cleanup
- **Organization**: Implement model creation in functions that return complete objects
- **Output**: Use specific output directories for generated files

## Project Structure
- `src/`: Contains runnable Python scripts
- `openmc/`: Contains OpenMC Python package
