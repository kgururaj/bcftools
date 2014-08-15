if [ "$#" -lt 2 ];
then
  echo "Needs at least 2 arguments  <reference_file>  '-l <file_list>'|'<vcf1> <vcf2> ...'";
  exit
fi
reference_file=$1;
shift
./bcftools merge  --merge-config=config.sample --gatk -i BaseQRankSum:median,MQ:median,MQ0:median,ClippingRankSum:median,MQRankSum:median,ReadPosRankSum:median,DP:sum -m all --gvcf --reference=${reference_file} $@
