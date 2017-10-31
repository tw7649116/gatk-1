package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

final class TempFineTuneMultipleAlignment {

    /**
     * Splits input alignments into ones that are intended to be used for chimeric alignments and
     * ones that are not reliable for that purpose,
     * based on provided MQ and unique ref/read span length threshold.
     */
    private static Tuple2<List<AlignmentInterval>, List<AlignmentInterval>> classifyAlignments(final Iterator<AlignmentInterval> iterator,
                                                                                               final int mapQualThreshold,
                                                                                               final int uniqueRefSpanThreshold,
                                                                                               final int uniqueReadSpanThreshold) {

        final List<AlignmentInterval> good = new ArrayList<>(10); // 10 is a blunt guess
        final List<AlignmentInterval> bad  = new ArrayList<>(10);

        AlignmentInterval current = iterator.next();
        while ( iterator.hasNext() ) {
            final AlignmentInterval next = iterator.next();
            final AlnPairUniqueLength alnPairUniqueLength = new AlnPairUniqueLength(current, next);

            if (alignmentIsNonInformative(current.mapQual, mapQualThreshold,
                    alnPairUniqueLength.oneUniqRefLen, uniqueRefSpanThreshold, alnPairUniqueLength.oneUniqReadLen, uniqueReadSpanThreshold)) {
                bad.add(current);
                current = next;
            } else if (alignmentIsNonInformative(next.mapQual, mapQualThreshold,
                    alnPairUniqueLength.twoUniqRefLen, uniqueRefSpanThreshold, alnPairUniqueLength.twoUniqReadLen, uniqueReadSpanThreshold)) {
                bad.add(next);
            } else {
                good.add(current);
                current = next;
            }
        }

        return new Tuple2<>(good, bad);
    }

    private static boolean alignmentIsNonInformative(final int mapQ, final int mapqThresholdInclusive,
                                                     final int uniqRefSpanLen, final int uniqueRefSpanThreshold,
                                                     final int uniqReadSpanLen, final int uniqueReadSpanThreshold) {
        return mapQ < mapqThresholdInclusive
                || uniqRefSpanLen < uniqueRefSpanThreshold
                || uniqReadSpanLen < uniqueReadSpanThreshold;

    }

    /**
     * For representing unique reference span sizes and read consumption length values of two neighboring
     * alignment intervals of a particular contig.
     * Fields are mostly useful, for now, for filtering alignments.
     */
    private static final class AlnPairUniqueLength {
        final int oneUniqRefLen;
        final int oneUniqReadLen;
        final int twoUniqRefLen;
        final int twoUniqReadLen;

        AlnPairUniqueLength(final AlignmentInterval one, final AlignmentInterval two) {
            Utils.validateArg(one.startInAssembledContig <= two.startInAssembledContig,
                    "assumption that input alignments are order along read is violated");

            final int overlapOnRefSpan = AlignmentInterval.overlapOnRefSpan(one, two);
            final int overlapOnRead = AlignmentInterval.overlapOnContig(one, two);

            if (overlapOnRead == 0) {
                oneUniqRefLen = one.referenceSpan.size() - overlapOnRefSpan;
                twoUniqRefLen = two.referenceSpan.size() - overlapOnRefSpan;
                oneUniqReadLen = one.endInAssembledContig - one.startInAssembledContig + 1;
                twoUniqReadLen = two.endInAssembledContig - two.startInAssembledContig + 1;
            } else {
                // TODO: 10/16/17 hardclip offset
                final int i = SvCigarUtils.computeAssociatedDistOnRef(one.cigarAlong5to3DirectionOfContig, two.startInAssembledContig, overlapOnRead);
                final int j = SvCigarUtils.computeAssociatedDistOnRef(two.cigarAlong5to3DirectionOfContig, two.startInAssembledContig, overlapOnRead);
                oneUniqRefLen = one.referenceSpan.size() - Math.max(i, overlapOnRefSpan);
                twoUniqRefLen = two.referenceSpan.size() - Math.max(j, overlapOnRefSpan);
                oneUniqReadLen = one.endInAssembledContig - one.startInAssembledContig + 1 - overlapOnRead;
                twoUniqReadLen = two.endInAssembledContig - two.startInAssembledContig + 1 - overlapOnRead;
            }
        }
    }
}
