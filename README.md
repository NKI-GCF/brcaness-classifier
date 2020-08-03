
# The docker image needs to be built to run thhe docker container. This takes a few hours, it is required once.

docker build --tag=ovabrca:latest .
# dependent on your configuration you may have to add the flag '--network=host' to the above.


# the docker container requires a co




data_dir= #/net/NGSarchive/RTAdump/191112_D00645_0382_ACE43GANXX_RUN993/Data/Intensities/BaseCalls/Project_RUN993_Simone_Koole_5610/
bwaindex_dir=/net/NGSanalysis/ref/Homo_sapiens.GRCh38/index/bwa/
path="`pwd`"



docker run --rm -t -i -u $UID:$GROUPS \
  -v $path/output/:/output:rw \
  -v $data_dir:/input:ro \
  -v $bwaindex_dir:/ref:ro \
  -v $path/config:/config:ro \
  -v $path/:/app:rw \
  ovabrca:0.6 bash

docker run --rm -u $UID:$GROUPS \
  -v $path/output/:/output:rw \
  -v $data_dir:/input:ro \
  -v $bwaindex_dir:/ref:ro \
  -v $path/config:/config:ro \
  ovabrca:0.6

