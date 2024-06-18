function import_config {
while read var value
do
if [ -z "$var" ]; then
	echo "mannaggia"
	else
	export "$var"="$value"
fi
done < $config_path
}

function check_defaults {
	if [ "$stats_path" == "default" ]; then
		stats_path="$bam_path/STATS"
	fi

	if [ "$motif_path" == "default" ]; then
		motif_path="$stats_path/MOTIF"
	fi

	if [ "$motif_count_path" == "default" ]; then
		motif_count_path="$stats_path/MOTIF_COUNTS"
	fi

	if [ "$readl_path" == "default" ]; then
		readl_path="$stats_path/READLENGTH_COUNTS"
	fi

	if [ "$raw_readl_path" == "default" ]; then
		raw_readl_path="$fastq_path/READL_RAW"
	fi
}
