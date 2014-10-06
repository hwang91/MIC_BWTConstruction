//


int main(int argc, char** argv){

	ParseArguments();
	LoadReadsSet();

	for(prefix = 0; prefix < 256; prefix++){
		AllocMemory();
		Partition();
		AmericanFlagSort();
		OutputBWT();
		FreeMemory();
	}

	return 0;
}