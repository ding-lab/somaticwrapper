# Be sure to do the following first:
# docker login cgc-images.sbgenomics.com

IMAGE="cgc-images.sbgenomics.com/m_wyczalkowski/somatic-wrapper:cwl"
#docker push mwyczalkowski/somatic-wrapper:cwl
docker push $IMAGE

