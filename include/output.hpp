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

 /*
 * @file output.hpp
 * Created: Nov 06, 2017 Mon
 * @brief Class for handling SMR output
 */

#pragma once

#include "report_fastx.h"
#include "report_fx_other.h"
#include "report_blast.h"
#include "report_denovo.h"
#include "report_sam.h"
#include "report_biom.h"

// forward
struct Index;
class References;
class Read;
class Refstats;
struct Readstats;
struct Runopts;
class KeyValueDatabase;
class Readfeed;

void writeReports(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts);

class Output {
public:
	ReportFastx fastx; // fastx aligned report
	ReportFxOther fx_other; // fastx non-aligned report
	ReportBlast blast;
	ReportDenovo denovo;
	ReportSam sam;
	ReportBiom biom;

	Output(Readfeed& readfeed, Runopts& opts);
	//~Output();

private:
	void init(Readfeed& readfeed, Runopts& opts);

}; // ~class Output