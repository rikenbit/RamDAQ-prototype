cd docker-cwltool-with-bash/
docker build -t cwltool-with-bash:1.0.20180809224403 .
docker tag cwltool-with-bash:1.0.20180809224403 yuifu/cwltool-with-bash:1.0.20180809224403
docker push yuifu/cwltool-with-bash:1.0.20180809224403
