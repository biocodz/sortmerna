/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent No�      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mika�l Salson    mikael.salson@lifl.fr
			   H�l�ne Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

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
	openfw(opts.dbg_level); // open output files for writing
	is_zip = (opts.zip_out == 1) || (readfeed.orig_files[0].isZip && opts.zip_out == -1);
	if (is_zip)	init_zip();
} // ~ReportFastx::init

void ReportFastx::append(const uint32_t& id, std::vector<Read>& reads, const Runopts& opts, const bool& is_last)
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
		for (uint32_t i = 0, idx = 0; i < reads.size(); ++i, idx = 0)
		{
			// 1 output file a (aligned reads)
			if (base.num_out == 1) {
				if (opts.is_paired_out) {
					if (reads[0].is_hit && reads[1].is_hit)
						idx = id;
					else
						continue; // to other 
				}
				else if (opts.is_paired_in || reads[i].is_hit)
					idx = id;
				else
					continue;
			}
			// 2 output files ap,as (sout) | af,ar (out2)
			else if (base.num_out == 2) {
				if (opts.is_out2) {
					if (opts.is_paired_out) {
						if (reads[0].is_hit && reads[1].is_hit) {
							idx = id * base.num_out + i;
						}
						else 
							break; // if not both aligned -> non-aligned
					}
					else if (opts.is_paired_in || reads[i].is_hit)
						idx = id * base.num_out + i;
					else
						continue; // ignore non-aligned
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
				base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx], is_last, opts.dbg_level);
			else
				base.write_a_read(fsv[idx], reads[i], opts.dbg_level);
		}
	}//~if paired
	// non-paired
	else
	{
		// the read was accepted - output
		if (reads[0].is_hit)
		{
			if (is_zip)
				base.write_a_read(fsv[id], reads[0], vstate_out[id], vzlib_out[id], is_last, opts.dbg_level);
			else
				base.write_a_read(fsv[id], reads[0], opts.dbg_level);
		}
	}
} // ~ReportFasta::append

ReportFxBase& ReportFastx::getBase()
{
	return base;
}