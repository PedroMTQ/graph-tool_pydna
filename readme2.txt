These modules implement the package Graph-tool.
Since there wasn't an improvement in efficiency it's probably better not to merge these modules into pydna.
Instead, a separate folder was created so that the code can be seen and evaluated.
This folder also contains several profiling files with both NetworkX and Graph-tool.
The .py files within this folder contain the tools necessary for profiling and timing of the code, so replication of the process should be trivial.
