package org.broadinstitute.hellbender.tools.walkers.mutect;

public enum OverlappingReadPairStrategy {
    KEEP_HIGHER_QUALITY_READ(true), KEEP_BOTH(false), ADD_LIKELIHOODS(true), ADD_PHRED_TEN(true), SET_TO_PHRED_FORTY_FIVE(true);

    private final boolean discardWorseRead;

    OverlappingReadPairStrategy(final boolean discardWorseRead) {
        this.discardWorseRead = discardWorseRead;
    }

    public boolean isDiscardWorseRead() { return discardWorseRead; }
}
