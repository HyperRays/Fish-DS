    sudo apt update

    sudo apt install -y apt-transport-https ca-certificates software-properties-common

    # add docker to known linux keyring
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | gpg --dearmor > docker-repo-keyring.gpg
    sudo mv ./docker-repo-keyring.gpg /etc/apt/trusted.gpg.d/
    sudo chown root:root /etc/apt/trusted.gpg.d/docker-repo-keyring.gpg
    sudo chmod ugo+r /etc/apt/trusted.gpg.d/docker-repo-keyring.gpg
    sudo chmod go-w /etc/apt/trusted.gpg.d/docker-repo-keyring.gpg
    sudo add-apt-repository -y "deb [arch=amd64] https://download.docker.com/linux/ubuntu focal stable"
    apt-cache policy docker-ce


    # install docker packages
    sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

    # Install newuidmap & newgidmap binaries
    apt install -y uidmap



    # rootless docker additional installation script
    # sudo apt-get install -y docker-ce-rootless-extras
    # dockerd-rootless-setuptool.sh install


    sh start_docker.sh