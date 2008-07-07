#!/usr/bin/perl

use Test::More tests => 108;
use Algorithm::ClusterPoints;

use Data::Dumper;

for my $dim (2, 3, 5) {
    for my $n (1, 5, 30, 200) {
        my @points = map rand, 1..$n*$dim;
        my $clp = Algorithm::ClusterPoints->new(radius => 1, ordered => 1, minimum_size => 1, dimension => $dim);
        $clp->add_points(@points);
        for my $ir (1, 10, 100) {
            my $r = 1/$ir;
            $clp->radius($r);
            for my $min_size (1, 2, 10) {
                my @clusters = $clp->clusters_ix;
                my @bfclusters = $clp->brute_force_clusters_ix;
                # print STDERR Data::Dumper->Dump([$r, \@clusters, \@bfclusters], [qw(r clusters bfclusters)]);
                is_deeply(\@clusters, \@bfclusters, "dim: $dim, n: $n, ir: $ir");
            }
        }
    }
}

