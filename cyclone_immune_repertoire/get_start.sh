json="cyclone_immune_repertoire_input.json"
java -jar /jdfssz2/ST_SUPERCELLS/P22Z10200N0739/liuyi/sorfware/wdl/cromwell-40.jar run \
	-i $json\
	/jdfssz2/ST_SUPERCELLS/P22Z10200N0739/liuyi/scirpt/cyclone_immune_repertoire/cyclone_immune_repertoire_v1.wdl
