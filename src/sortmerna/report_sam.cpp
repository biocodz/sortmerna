#include "report_sam.h"
#include "common.hpp"
#include "options.hpp"
#include "read.hpp"
#include "references.hpp"
#include "refstats.hpp"
#include "readfeed.hpp"

ReportSam::ReportSam(Runopts& opts) : Report(opts) {}

ReportSam::ReportSam(Readfeed& readfeed, Runopts& opts) : ReportSam(opts)
{
	init(readfeed, opts);
}

void ReportSam::init(Readfeed& readfeed, Runopts& opts)
{
	fv.resize(readfeed.num_splits);
	fsv.resize(readfeed.num_splits);
	is_zip = readfeed.orig_files[0].isZip;
	// WORKDIR/out/aligned_0_PID.sam
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + pid_str : "";
		std::string gz = is_zip ? ".gz" : "";
		fv[i] = opts.aligned_pfx.string() + sfx1 + sfx2 + ext + gz;
		fv[i] = opts.aligned_pfx.string() + "_denovo" + sfx1 + sfx2 + ext + gz;
		openfw(i);
	}
	if (is_zip) init_zip();
}

void ReportSam::append(int id, Read& read, References& refs, Runopts& opts)
{
	std::stringstream ss;
	if (read.is03) read.flip34();

	// read did not align, output null string
	if (opts.is_print_all_reads && read.alignment.alignv.size() == 0)
	{
		// (1) Query
		ss << read.getSeqId();
		ss << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		return;
	}

	// read aligned, output full alignment
	// iterate read alignments
	for (int i = 0; i < read.alignment.alignv.size(); ++i)
	{
		if (read.alignment.alignv[i].index_num == refs.num
			&& read.alignment.alignv[i].part == refs.part)
		{
			// (1) Query
			ss << read.getSeqId();
			// (2) flag Forward/Reversed
			if (!read.alignment.alignv[i].strand) ss << "\t16\t";
			else ss << "\t0\t";
			// (3) Subject
			ss << refs.buffer[read.alignment.alignv[i].ref_num].id;
			// (4) Ref start
			ss << "\t" << read.alignment.alignv[i].ref_begin1 + 1;
			// (5) mapq
			ss << "\t" << 255 << "\t";
			// (6) CIGAR
			// output the masked region at beginning of alignment
			if (read.alignment.alignv[i].read_begin1 != 0)
				ss << read.alignment.alignv[i].read_begin1 << "S";

			for (int c = 0; c < read.alignment.alignv[i].cigar.size(); ++c)
			{
				uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
				uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
				ss << length;
				if (letter == 0) ss << "M";
				else if (letter == 1) ss << "I";
				else ss << "D";
			}

			auto end_mask = read.sequence.size() - read.alignment.alignv[i].read_end1 - 1;
			// output the masked region at end of alignment
			if (end_mask > 0) ss << end_mask << "S";
			// (7) RNEXT, (8) PNEXT, (9) TLEN
			ss << "\t*\t0\t0\t";
			// (10) SEQ

			if (read.alignment.alignv[i].strand == read.reversed) // XNOR
				read.revIntStr();
			ss << read.get04alphaSeq();
			// (11) QUAL
			ss << "\t";
			// reverse-complement strand
			if (read.quality.size() > 0 && !read.alignment.alignv[i].strand)
			{
				std::reverse(read.quality.begin(), read.quality.end());
				ss << read.quality;
			}
			else if (read.quality.size() > 0) // forward strand
			{
				ss << read.quality;
				// FASTA read
			}
			else ss << "*";

			// (12) OPTIONAL FIELD: SW alignment score generated by aligner
			ss << "\tAS:i:" << read.alignment.alignv[i].score1;
			// (13) OPTIONAL FIELD: edit distance to the reference
			uint32_t mismatches = 0;
			uint32_t gaps = 0;
			uint32_t mid = 0;
			read.calcMismatchGapId(refs, i, mismatches, gaps, mid);
			ss << "\tNM:i:" << mismatches + gaps << "\n";
		}
	} // ~for read.alignments
	if (is_zip) {
		auto ret = vzlib_out[id].defstr(ss.str(), fsv[id]); // Z_STREAM_END | Z_OK - ok
		if (ret < Z_OK || ret > Z_STREAM_END) {
			ERR("Failed deflating readstring: ", ss.str(), " zlib status: ", ret);
		}
	}
	else {
		fsv[id] << ss.str();
	}
} // ~ReportSam::append

void ReportSam::write_header(Runopts& opts)
{
	std::stringstream ss;
	ss << "@HD\tVN:1.0\tSO:unsorted\n";

	// TODO: this line is taken from "Index::load_stats". To be finished (20171215).
#if 0
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); index_num++)
	{
		//@SQ header
		if (opts.yes_SQ) acceptedsam << "@SQ\tSN:" << s << "\tLN:" << len_seq << "\n";
		// number of nucleotide sequences in the reference file
		uint32_t num_sq = 0;
		stats.read(reinterpret_cast<char*>(&num_sq), sizeof(uint32_t));

		// loop through each @SQ
		for (uint32_t j = 0; j < num_sq; j++)
		{
			// length of the sequence id
			uint32_t len_id = 0;
			stats.read(reinterpret_cast<char*>(&len_id), sizeof(uint32_t));
			// the sequence id string
			std::string s(len_id + 1, 0); // AK
			std::vector<char> vs(s.begin(), s.end());
			stats.read(reinterpret_cast<char*>(&vs[0]), sizeof(char) * len_id);
			// the length of the sequence itself
			uint32_t len_seq = 0;
			stats.read(reinterpret_cast<char*>(&len_seq), sizeof(uint32_t));
			//		 @SQ header
			if (opts.yes_SQ) acceptedsam << "@SQ\tSN:" << s << "\tLN:" << len_seq << "\n";
		} // ~for
	} // ~for
#endif
	ss << "@PG\tID:sortmerna\tVN:1.0\tCL:" << opts.cmdline << std::endl;
	if (is_zip) {
		auto ret = vzlib_out[0].defstr(ss.str(), fsv[0]); // Z_STREAM_END | Z_OK - ok
		if (ret < Z_OK || ret > Z_STREAM_END) {
			ERR("Failed deflating readstring: ", ss.str(), " zlib status: ", ret);
		}
	}
	else
		fsv[0] << ss.str();
} // ~ReportSam::write_header