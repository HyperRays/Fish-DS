FROM intel/oneapi-hpckit

RUN apt install -y make
WORKDIR /app/
ENTRYPOINT ["/bin/bash"]
