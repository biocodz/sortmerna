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

#pragma once

#include "report.h"
#include "report_fx_base.h"

// forward
class Read;

class ReportFastx: public Report
{
public:
	ReportFastx(Runopts& opts);
	ReportFastx(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	/*
    * called on each read or each pair of reads (if paired)
    * writes to aliged.fasta
    *
	* @param id  processing thread enumerator 0,1,2..
    * @param reads: 1 or 2 (paired) reads
	* @param opts
	* @param is_last   flags the passed read is the last
    */
	void append(const uint32_t& id, std::vector<Read>& reads, const Runopts& opts, const bool& is_last=false);
	ReportFxBase& getBase();

private:
	ReportFxBase base;
};
