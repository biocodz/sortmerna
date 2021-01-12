#include <filesystem>

#include "common.hpp"
#include "report_fastx.h"
#include "readfeed.hpp"
#include "options.hpp"
#include "read.hpp"

ReportFastx::ReportFastx(Runopts& opts): Report(opts), base() {}
ReportFastx::ReportFastx(Readfeed& readfeed, Runopts& opts): ReportFastx(opts)
{ 
	init(readfeed, opts); 
}

void ReportFastx::init(Readfeed& readfeed, Runopts& opts)
{
	base.init(opts);
	base.init(readfeed, opts, fv, fsv, opts.aligned_pfx.string(), pid_str);
	openfw(); // open output files for writing
	is_zip = readfeed.orig_files[0].isZip;
	if (is_zip)	init_zip();
} // ~ReportFastx::init

void ReportFastx::append(int id, std::vector<Read>& reads, Runopts& opts, bool is_last)
{
	if (opts.is_paired) {
		// validate the reads are indeed paired in case of two reads files
		if (opts.readfiles.size() == 2
			&& (reads[0].read_num != reads[1].read_num
				|| reads[0].readfile_idx == reads[1].readfile_idx))
		{
			ERR("Paired validation failed: reads[0].id= ", reads[0].id, " reads[0].read_num = ",
				reads[0].read_num, " reads[0].readfile_idx= ", reads[0].readfile_idx,
				" reads[1].id=", reads[1].id, " reads[1].read_num = ", reads[1].read_num,
				" reads[1].readfile_idx = ", reads[1].readfile_idx);
			exit(EXIT_FAILURE);
		}

		if (!reads[0].is_hit && !reads[1].is_hit)
			return; // neiher is aligned - nothing to write

		// caclulate the index of the output file to write to
		for (int i = 0, idx = 0; i < reads.size(); ++i)
		{
			// 1 output file a (aligned reads)
			if (base.num_out == 1) {
				if (reads[i].is_hit || opts.is_paired_in)
					idx = id;
				else
					continue;
			}
			// 2 output files ap,as (sout) | af,ar (out2)
			else if (base.num_out == 2) {
				if (opts.is_out2) {
					if (opts.is_paired_out && !(reads[0].is_hit && reads[1].is_hit))
						break; // if not both aligned -> non-aligned
					else if (opts.is_paired_in || reads[i].is_hit)
						idx = id * base.num_out + i;
				}
				else if (opts.is_sout) {
					if (reads[0].is_hit && reads[1].is_hit)
						idx = id * base.num_out; // both to 'ap' [0]
					else if (reads[i].is_hit)
						idx = id * base.num_out + 1; // hit to 'as' [1]
					else
						continue; // ignore non-aligned read
				}
			}
			// 4 output files: apf, apr, asf, asr
			// being here means both is_out2 & is_sout were set => No paired_in/out
			else if (base.num_out == 4) {
				if (reads[0].is_hit && reads[1].is_hit)
					idx = id * base.num_out + i; // 0 -> apf, 1 -> apr
				else if (reads[i].is_hit)
					idx = id * base.num_out + i + 2; // 0 -> asf, 1 -> asr
				else
					continue; // ignore a non-aligned singleton
			}
			else {
				ERR("min number of output files is 1, max number of output files is 4. The current value is ", base.num_out);
				exit(1);
			}

			if (is_zip)
				base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx], is_last);
			else
				base.write_a_read(fsv[idx], reads[i]);
		}
	}//~if paired
	// non-paired
	else
	{
		// the read was accepted - output
		if (reads[0].is_hit)
		{
			if (is_zip)
				base.write_a_read(fsv[0], reads[0], vstate_out[0], vzlib_out[0], is_last);
			else
				base.write_a_read(fsv[0], reads[0]);
		}
	}
} // ~ReportFasta::append

void ReportFastx::merge(int num_splits)
{
	for (int i = 0; i < base.num_out; ++i) {
		openfw(i);
		for (int j = 1; j < num_splits; ++j) {
			auto idx = i + j * base.num_out;
			openfr(idx);
			fsv[i] << fsv[idx].rdbuf();
			INFO("merged ", fv[idx], " -> ", fv[i]);
			closef(idx);
			std::filesystem::remove(fv[idx]);
			INFO("deleted ", fv[idx]);
		}
		closef(i);
		strip_path_sfx(fv[i]);
	}
}

ReportFxBase& ReportFastx::getBase()
{
	return base;
}