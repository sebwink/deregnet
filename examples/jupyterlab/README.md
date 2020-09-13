# Running DeRegNet interactively withing Jupyter Lab

Depending on whether you use a named-user or floating license for the Gurobi solver
you need to run the run script from the *../../docker/named-user* or *../../docker/token-server*
respectively. In the following examples adapt the environment variable *GUROBI_LICENSE_MODE* variable
accordingly.

From this directory run:

```sh
export GUROBI_LICENSE_MODE=named-user
export DEREGNET_HOME=../..

$DEREGNET_HOME/docker/$GUROBI_LICENSE_MODE/run sebwink/deregnet:0.99.999 jupyter lab
```

Navigating to the displayed URL in your browser you can then execute the example notebook
*deregnet-example.ipynb*.

Your current working directory gets mounted into the Docker container and results of the benchmark 
script correspondingly will end up there after it finishes running.

For more control concerning Docker volumes and exposed ports, etc. you can have a look at the implementation of the
*docker/<license-mode>/run* scripts for [named-user](https://github.com/sebwink/deregnet/tree/master/docker/named-user) mode 
and [token server/floating](https://github.com/sebwink/deregnet/tree/master/docker/token-server) mode.
