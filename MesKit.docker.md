### Use MesKit with Docker Container 

Container technology makes applications more accessible and IDEA is also dockerized and provides an image on [Docker Hub](https://cloud.docker.com/repository/docker/niinleslie/meskit/builds).


#### Prerequisites


In order to run MesKit with docker containers, you'll need docker installed.

* [Windows](https://docs.docker.com/windows/started)
* [OS X](https://docs.docker.com/mac/started/)
* [Linux](https://docs.docker.com/linux/started/)

#### Usage

To get this MesKit image from Docker Hub,  you can execute the following command in shell :

```shell
$ docker pull niinleslie/meskit
```

After pulling MesKit image checked by `docker images`, you should run the `docker run` command below to initiate the docker container and make sure that the port is suitable for shiny-server's configuration: 

```shell
$ docker run -p 3838:3838 niinleslie/MesKit shiny-server
```

If the command was workable, the docker container would start running and the shiny server of MesKit would be deployed. To , you could see the following messages:

```shell
$ docker run -p 3838:3838 niinleslie/mesKit shiny-server
...
[2019-02-20T02:04:56.635] [INFO] shiny-server - Starting listener on http://<IP_address>:3838
...
```

And the third `[INFO]`message showed above gives a local link (` http://localhost:3838` or `http://<IP_address>:3838`)to use MesKit with port 3838. Simply copy or open the link by your web browser and if the platform is successfully loaded, that means you can agreeably enjoy your MesKit trip. 

Finally, if you have got all your work done, you can close your web browser and turn off the server simply by `Ctrl+C`+ `Enter`. And now you can `exit` the docker container and continue your analysis with downloaded data files from MesKit.
