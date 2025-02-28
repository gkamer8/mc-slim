
# Setup

First, build the Docker container:

```
docker build -t my_openmc_custom .
```

Now run it with the appropriate volume mount so you can easily edit the source code:

```
docker run -it --name=my_openmc_custom -v ./src:/root/src my_openmc_custom
```

If you exit, you can reenter with:

```
docker start -i my_openmc_custom
```

Inside the docker container, you'll be placed in `~`, where you can enter `src` and run relevant python files, like `python3 pincell.py`.

## VSCode

The openmc package referenced in the code exists only in the docker container - not necessarily on a host machine. Therefore, to give VSCode linting in some dev setups, it is useful to copy the OpenMC module from [GitHub](https://github.com/openmc-dev/openmc) into the repo (in particular, the OpenMC Python package with `__init__.py` at its root located in `/openmc` in the repo).
