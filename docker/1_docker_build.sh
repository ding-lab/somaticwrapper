# Be sure to do the following first:
# docker login cgc-images.sbgenomics.com

#IMAGE="mwyczalkowski/somatic-wrapper:cwl"
IMAGE="cgc-images.sbgenomics.com/m_wyczalkowski/somatic-wrapper:cwl-dev"

docker build -t $IMAGE -f Dockerfile .
