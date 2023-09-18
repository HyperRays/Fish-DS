# Fish-DS


# Installation

Make sure you have docker installed. If not here is the link: \
https://www.docker.com/

Run `docker compose up` this will set up the docker container and start it. \
Since the image is quite large it may take up to 20 min to complete.

If you have the vscode docker extension installed you can attach a shell and run `bash`

If not, go to your cli and find the container using `docker ps` \
note down the container id then run the command `docker exec -it bash <container-id>`

Using `cd`, navigate to `fish-basic`. To compile the program, run `make` \
To run the program execute `mpirun ./fish`

For further infos read the document written by Mr. Liebendörfer `fish-basic/readme.md`


