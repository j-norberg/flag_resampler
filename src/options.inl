
struct Options
{
	// bool test = false;
	bool show_version = false;
	bool show_usage = false;

	std::string in_file;
	std::string out_file;
	int quality = 100;
	Writer::OutFormat _format = Writer::eFmtFloat;

	int out_sr = 44100;
};

const static char* next_arg(int argc, const char** argv, int& index, const char* error_msg)
{
	if (index < (argc - 1))
	{
		++index;
		return argv[index];
	}

	puts(error_msg);
	return nullptr;
}

bool handle_flag(Options& o, int argc, const char** argv, int& index)
{
	std::string arg(argv[index]);

#if 0
	if (arg == "-t" || arg == "--test")
	{
		o.test = true;
		return true;
	}
#endif

	if (arg == "-v" || arg == "--version")
	{
		o.show_version = true;
		return true;
	}

	if (arg == "-i" || arg == "--input")
	{
		o.in_file = next_arg(argc, argv, index, "error: -i should be followed by filename\n");
		return o.in_file.empty() == false;
	}

	if (arg == "-o" || arg == "--output")
	{
		o.out_file = next_arg(argc, argv, index, "error: -o should be followed by filename\n");
		return o.out_file.empty() == false;
	}

	if (arg == "-r" || arg == "--sample_rate")
	{
		const char* sr_str = next_arg(argc, argv, index, "error: -r should be followed by sample rate\n");
		o.out_sr = sr_str ? atoi(sr_str) : 0; // set to 0 to indicate issue
		return o.out_sr > 0;
	}

	if (arg == "-q" || arg == "--quality")
	{
		const char* sr_str = next_arg(argc, argv, index, "error: -q should be followed by quality percentage\n");
		o.quality = sr_str ? atoi(sr_str) : 0; // set to 0 to indicate issue
		return o.quality > 0;
	}

	if (arg == "-f" || arg == "--format")
	{
		const char* format_error_message = "error: -f should be followed by format (16, 24, f32, f64)\n";
		const char* fmt_str = next_arg(argc, argv, index, format_error_message);
		if (fmt_str == nullptr)
			return false;

		if (0 == strcmp(fmt_str, "16"))
		{
			o._format = Writer::eFmtInt16Dithered;
			return true;
		}

		if (0 == strcmp(fmt_str, "24"))
		{
			o._format = Writer::eFmtInt24Dithered;
			return true;
		}

		if (0 == strcmp(fmt_str, "f32"))
		{
			o._format = Writer::eFmtFloat;
			return true;
		}

		if (0 == strcmp(fmt_str, "f64"))
		{
			o._format = Writer::eFmtDouble;
			return true;
		}

		puts(format_error_message);
		return false;
	}

	printf("unhandled flag %s\n\n", arg.c_str());
	return false;
}


Options get_options(int argc, const char** argv)
{
	// fixme more reporting about invalid flags
	Options o;

	for (int i = 1; i < argc; ++i)
		if (!handle_flag(o, argc, argv, i))
		{
			o.out_file.clear(); // make invalid
			break;
		}

	// when do we need to show usage?
//	if (!o.test)
	{
		if (o.in_file.empty() || o.out_file.empty())
			o.show_usage = true;
	}

	return o;
};

