# Fish-DS


# Installation

Make sure you have docker installed. 

run `docker compose up` this will set up the docker container and start it
since the image is quite large it may take up to 20 min

If you have the vscode docker extension installed you can attach a shell and run `bash`

If not, go to your cli and find the container using `docker ps`
note down the container id then run the command `docker exec -it bash <container-id>`

using `cd` navigate to `fish-basic` and to compile the program run `make`
to run the program run `mpirun ./fish`

for further infos read the document written by Mr. Liebendörfer `fish-basic/readme.md`


