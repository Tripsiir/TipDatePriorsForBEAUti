for i in ./templates/*
do
	name=${i##*/}

	if [ ${name:0:2} -eq '01' -o ${name:0:2} -eq '02' ]
	then
		continue  ### resumes iteration of an enclosing for loop ###
	fi

	if grep -q realspace <<<${name}; then
		python ../sample-tip-dates-in-BEAUti-xml.py ./templates/02-tip-dates-no-prior.xml --sequences restrictedSequenceList.txt --output validationOutput/updated-log-normal-realspace.xml -d log-normal --realspace
	    continue
	fi

	if grep -q '1x' <<<${name}; then
	    python ../sample-tip-dates-in-BEAUti-xml.py ./templates/02-tip-dates-no-prior.xml --sequences restrictedSequenceList.txt --output validationOutput/updated-${name} -d 1/x
	    continue
	fi

	dist=$(echo ${name} | sed -r 's/[0-9]{2}-([[:alnum:]-]*).xml/\1/g')
	python ../sample-tip-dates-in-BEAUti-xml.py ./templates/02-tip-dates-no-prior.xml --sequences restrictedSequenceList.txt --output validationOutput/updated-${name} -d ${dist}

done
