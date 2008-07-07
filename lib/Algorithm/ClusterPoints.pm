package Algorithm::ClusterPoints;

our $VERSION = '0.04';

use strict;
use warnings;

use constant sqr2 => sqrt(2);
use constant isqr2 => 1/sqr2;

use POSIX qw(floor ceil);
use List::Util qw(max min);
use Carp;

use Data::Dumper;

my $packing = ($] >= 5.008 ? 'w*' : 'V*');

sub new {
    @_ & 1 or croak 'Usage: Algorithm::ClusterPoints->new(%options)';
    my ($class, %opts) = @_;

    my $dimension = delete $opts{dimension};
    $dimension = 2 unless defined $dimension;
    $dimension < 1 and croak "positive dimension required";
    my $radius = delete $opts{radius};
    my $minimum_size = delete $opts{minimum_size};
    my $ordered = delete $opts{ordered};

    %opts and croak "unknown constructor option(s) '".join("', '", sort keys %opts)."'";

    my $self = bless { radius => 1.0,
                       minimum_size => 1,
                       ordered => 0,
                       dimension => $dimension,
                       coords => [ map [], 1..$dimension ],
                     }, $class;

    $self->radius($radius) if defined $radius;
    $self->minimum_size($minimum_size) if defined $minimum_size;
    $self->ordered($ordered) if defined $ordered;
    $self;
}

sub add_point {
    my $self = shift;
    my $dimension = $self->{dimension};
    @_ % $dimension and croak 'Usage: $clp->add_point($x0, $y0, ... [,$x1, $y1 ...])';
    delete $self->{_clusters};
    my $ix = @{$self->{coords}[0]};
    while (@_) {
        push @$_, shift
            for (@{$self->{coords}});
    }
    $ix;
}

*add_points = \&add_point;

sub point_coords {
    @_ == 2 or croak 'Usage: $clp->point_coords($index)';
    my ($self, $ix) = @_;
    my $top = $#{$self->{coords}[0]};
    croak "point index $ix out of range [0, $top]"
        if ($ix > $top or $ix < 0);
    return $self->{coords}[0][$ix]
        if $self->{dimension} == 1;
    wantarray or croak 'method requires list context';
    map $_->[$ix], @{$self->{coords}};
}

sub reset { delete shift->{_clusters} }

sub radius {
    @_ > 2 and croak 'Usage: $clp->radius([$new_radius])';
    my $self = shift;
    if (@_) {
        my $radius = shift;
        $radius > 0.0 or croak 'positive radius required';
        $self->{radius} = $radius;
        delete $self->{_clusters};
    }
    $self->{radius};
}

sub ordered {
    @_ > 2 and croak 'Usage: $clp->ordered([$ordered])';
    my $self = shift;
    if (@_) {
        $self->{ordered} = !!shift;
        delete $self->{_clusters};
    }
    $self->{ordered};
}

sub minimum_size {
    @_ > 2 and croak 'Usage: $clp->minimum_size([$size])';
    my $self = shift;
    if (@_) {
        my $minimum_size = shift;
        $minimum_size > 0 or croak 'positive minimum_size required';
        $self->{minimum_size} = $minimum_size;
        delete $self->{_clusters};
    }
    $self->{minimum_size};
}

sub clusters {
    my $self = shift;
    my $clusters = $self->{_clusters} ||= $self->_make_clusters_ix;
    my $ax = $self->{x};
    my $ay = $self->{y};
    my $coords = $self->{coords};
    my @out;
    for my $cluster (@$clusters) {
        my @cluster_coords;
        for my $ix (@$cluster) {
            push @cluster_coords, map $_->[$ix], @$coords;
        }
        push @out, \@cluster_coords;
    }
    @out;
}

sub clusters_ix {
    my $self = shift;
    my $clusters = $self->{_clusters} ||= $self->_make_clusters_ix;
    # dup arrays:
    map [@$_], @$clusters;
}

sub _touch_2 {
    my ($radius, $c1, $c2, $ax, $ay) = @_;

    my $c1_xmin = min @{$ax}[@$c1];
    my $c2_xmax = max @{$ax}[@$c2];
    return 0 if $c1_xmin - $c2_xmax > $radius;

    my $c1_xmax = max @{$ax}[@$c1];
    my $c2_xmin = min @{$ax}[@$c2];
    return 0 if $c2_xmin - $c1_xmax > $radius;

    my $c1_ymin = min @{$ay}[@$c1];
    my $c2_ymax = max @{$ay}[@$c2];
    return 0 if $c1_ymin - $c2_ymax > $radius;

    my $c1_ymax = max @{$ay}[@$c1];
    my $c2_ymin = min @{$ay}[@$c2];
    return 0 if $c2_ymin - $c1_ymax > $radius;

    my $r2 = $radius * $radius;
    for my $i (@$c1) {
        for my $j (@$c2) {
            my $dx = $ax->[$i] - $ax->[$j];
            my $dy = $ay->[$i] - $ay->[$j];
            return 1 if ($dx * $dx + $dy * $dy <= $r2);
        }
    }
    0;
}

sub _touch {
    my ($radius, $c1, $c2, $coords) = @_;
    # print STDERR "touch($c1->[0], $c2->[0])\n";

    for my $coord (@$coords) {
        my $c1_min = min @{$coord}[@$c1];
        my $c2_max = max @{$coord}[@$c2];
        return 0 if $c1_min - $c2_max > $radius;

        my $c1_max = max @{$coord}[@$c1];
        my $c2_min = min @{$coord}[@$c2];
        return 0 if $c2_min - $c1_max > $radius;
    }

    my $r2 = $radius * $radius;
    for my $i (@$c1) {
        for my $j (@$c2) {
            my $sum = 0;
            for (@$coords) {
                my $delta = $_->[$i] - $_->[$j];
                $sum += $delta * $delta;
            }
            return 1 if $sum <= $r2;
        }
    }
    # print STDERR "don't touch\n";
    0;
}

sub _make_clusters_ix {
    my $self = shift;
    # print STDERR Dumper $self;
    $self->{dimension} == 2
        ? $self->_make_clusters_ix_2
        : $self->_make_clusters_ix_any;
}

sub _make_clusters_ix_2 {
    my $self = shift;

    my $radius = $self->{radius};
    my $scale = 1.00001*sqr2/$radius;

    $self->{dimension} == 2
        or croak 'internal error: _make_clusters_ix_2 called but dimension is not 2';

    my $ax = $self->{coords}[0];
    my $ay = $self->{coords}[1];
    @$ax or croak "points have not been added";

    my $xmin = min @$ax;
    my $xmax = max @$ax;
    my $ymin = min @$ay;
    my $ymax = max @$ay;

    my @fx = map { floor($scale * ($_ - $xmin)) } @$ax;
    my @fy = map { floor($scale * ($_ - $ymin)) } @$ay;

    my (%ifx, %ify, $c);
    $c = 1; $ifx{$_} ||= $c++ for @fx;
    $c = 1; $ify{$_} ||= $c++ for @fy;
    my %rifx = reverse %ifx;
    my %rify = reverse %ify;

    # print STDERR "radius: $radius\n";

    my %cell;
    # my %cellid;
    my $cellid = 1;
    for my $i (0..$#$ax) {
        my $cell = pack $packing => $ifx{$fx[$i]}, $ify{$fy[$i]};
        push @{$cell{$cell}}, $i;
        # $cellid{$cell} ||= $cellid++;
        # print STDERR "i: $i, x: $ax->[$i], y: $ay->[$i], fx: $fx[$i], fy: $fy[$i], ifx: $ifx{$fx[$i]}, ify: $ify{$fy[$i]}, cellid: $cellid{$cell}\n";
    }

    my %cell2cluster; # n to 1 relation
    my %cluster2cell; # when $cluster2cell{$foo} does not exists

    while(defined (my $cell = each %cell)) {
        my %cluster;
        my ($ifx, $ify) = unpack $packing => $cell;
        my $fx = $rifx{$ifx};
        my $fy = $rify{$ify};
        for my $dx (-2, -1, 0, 1, 2) {
            my $ifx = $ifx{$fx + $dx};
            defined $ifx or next;
            my $filter = 6 - $dx * $dx;
            for my $dy (-2, -1, 0, 1, 2) {
                # next if $dx * $dx + $dy * $dy > 5;
                next if $dy * $dy > $filter;
                my $ify = $ify{$fy + $dy};
                defined $ify or next;
                my $neighbor = pack $packing => $ifx, $ify;
                my $cluster = $cell2cluster{$neighbor};
                if ( defined $cluster and
                     !$cluster{$cluster} and
                     _touch_2($radius, $cell{$cell}, $cell{$neighbor}, $ax, $ay) ) {
                    $cluster{$cluster} = 1;
                }
            }
        }
        if (%cluster) {
            my ($to, $to_cells);
            if (keys %cluster > 1) {
                my $max = 0;
                for (keys %cluster) {
                    my $cells = $cluster2cell{$_};
                    my $n = defined($cells) ? @$cells : 1;
                    if ($n > $max) {
                        $max = $n;
                        $to = $_;
                    }
                }
                delete $cluster{$to};
                $to_cells = ($cluster2cell{$to} ||= [$to]);
                for (keys %cluster) {
                    my $neighbors = delete $cluster2cell{$_};
                    if (defined $neighbors) {
                        push @$to_cells, @$neighbors;
                        $cell2cluster{$_} = $to for @$neighbors;
                    }
                    else {
                        push @$to_cells, $_;
                        $cell2cluster{$_} = $to;
                    }
                }
            }
            else {
                $to = each %cluster;
                $to_cells = ($cluster2cell{$to} ||= [$to]);
            }
            push @$to_cells, $cell;
            $cell2cluster{$cell} = $to;
        }
        else {
            $cell2cluster{$cell} = $cell;
        }
    }

    my @clusters;
    while (my ($cluster, $cells) = each %cluster2cell) {
        my @points = map @{delete $cell{$_}}, @$cells;
        if (@points >= $self->{minimum_size}) {
            @points = sort { $a <=> $b } @points if $self->{ordered};
            push @clusters, \@points;
        }
    }
    push @clusters, grep { @$_ >= $self->{minimum_size} } values %cell;

    @clusters = sort { $a->[0] <=> $b->[0] } @clusters if $self->{ordered};

    return \@clusters;
}


sub _delta_sphere {
    my $dimension = shift;
    my @delta_sphere = [];
    my $delta_top = ceil(sqrt($dimension));
    my @deltas = (-$delta_top..$delta_top);
    # print STDERR "deltas: @deltas\n";
    for (1..$dimension) {
        my @next;
        for my $ds (@delta_sphere) {
            push @next, map [@$ds, $_], @deltas;
        }
        @delta_sphere = @next;
    }
    my $delta_top2 = $delta_top * $delta_top;
    grep {
        my $sum = 0;
        for (@$_) {
            my $min = ($_ ? abs($_) - 1 : 0);
            $sum += $min * $min;
        }
        $sum < $delta_top2;
    } @delta_sphere;
}

my %delta_sphere; # cache
sub _make_clusters_ix_any {
    my $self = shift;

    my $radius = $self->{radius};
    my $dimension = $self->{dimension};
    my $dimension_top = $dimension - 1;
    my $scale = 1.00001*sqrt($dimension)/$radius;
    my $coords = $self->{coords};
    my $top = $#{$coords->[0]};
    $top >= 0 or croak "points have not been added";

    my (@fls, @ifls, @rifls);
    for my $coord (@$coords) {
        my $min = min @$coord;
        my @fl = map floor($scale * ($_ - $min)), @$coord;
        push @fls, \@fl;
        my %ifl;
        my $c = 1;
        $ifl{$_} ||= $c++ for @fl;
        push @ifls, \%ifl;
        my %rifl = reverse %ifl;
        push @rifls, \%rifl;
    }

    my %cell;
    for my $i (0..$top) {
        my $cell = pack $packing => map $ifls[$_]{$fls[$_][$i]}, 0..$dimension_top;
        push @{$cell{$cell}}, $i;
    }
    # print STDERR "\%cell:\n", Dumper [values %cell];

    my %cell2cluster; # n to 1 relation
    my %cluster2cell; # when $cluster2cell{$foo} does not exists

    my @delta_sphere = @{$delta_sphere{$dimension} ||= [_delta_sphere($dimension)]};
    # print STDERR "delta_sphere\n", Dumper \@delta_sphere;

    while(defined (my $cell = each %cell)) {
        my %cluster;
        my @ifl = unpack $packing => $cell;
        my @fl = map $rifls[$_]{$ifl[$_]}, 0..$dimension_top;

        for my $delta (@delta_sphere) {
            my @ifl = map { $ifls[$_]{$fl[$_] + $delta->[$_]} || next } 0..$dimension_top;
            # next if grep !defined, @ifl;
            my $neighbor = pack $packing => @ifl;
            my $cluster = $cell2cluster{$neighbor};
            if ( defined $cluster and
                 !$cluster{$cluster} and
                 _touch($radius, $cell{$cell}, $cell{$neighbor}, $coords) ) {
                $cluster{$cluster} = 1;
            }
        }
        if (%cluster) {
            my ($to, $to_cells);
            if (keys %cluster > 1) {
                my $max = 0;
                for (keys %cluster) {
                    my $cells = $cluster2cell{$_};
                    my $n = defined($cells) ? @$cells : 1;
                    if ($n > $max) {
                        $max = $n;
                        $to = $_;
                    }
                }
                delete $cluster{$to};
                $to_cells = ($cluster2cell{$to} ||= [$to]);
                for (keys %cluster) {
                    my $neighbors = delete $cluster2cell{$_};
                    if (defined $neighbors) {
                        push @$to_cells, @$neighbors;
                        $cell2cluster{$_} = $to for @$neighbors;
                    }
                    else {
                        push @$to_cells, $_;
                        $cell2cluster{$_} = $to;
                    }
                }
            }
            else {
                $to = each %cluster;
                $to_cells = ($cluster2cell{$to} ||= [$to]);
            }
            push @$to_cells, $cell;
            $cell2cluster{$cell} = $to;
        }
        else {
            $cell2cluster{$cell} = $cell;
        }
    }

    my @clusters;
    while (my ($cluster, $cells) = each %cluster2cell) {
        my @points = map @{delete $cell{$_}}, @$cells;
        if (@points >= $self->{minimum_size}) {
            @points = sort { $a <=> $b } @points if $self->{ordered};
            push @clusters, \@points;
        }
    }
    push @clusters, grep { @$_ >= $self->{minimum_size} } values %cell;

    @clusters = sort { $a->[0] <=> $b->[0] } @clusters if $self->{ordered};

    return \@clusters;
}

sub brute_force_clusters_ix {
    my $self = shift;
    $self->{dimension} == 2
        ? $self->brute_force_clusters_ix_2
        : $self->brute_force_clusters_ix_any
}

sub brute_force_clusters_ix_2 {
    @_ == 1 or croak 'Usage: $clp->brute_force_clusters_ix_2';
    my $self = shift;

    my $radius = $self->{radius};
    my $radius2 = $radius * $radius;

    $self->{dimension} == 2
        or croak "brute_force_clusters_ix_2 called but dimension is not equal to 2";
    
    my $ax = $self->{coords}[0];
    my $ay = $self->{coords}[1];
    @$ax or croak "points have not been added";

    my @ix = 0..$#$ax;
    my @clusters;
    while (@ix) {
        # print STDERR "\@ix 1: ".join("-", @ix).".\n";
        my @current = shift @ix;
        my @done;
        while (@ix) {
            # print STDERR "\@ix 2: ".join("-", @ix).".\n";
            my $ix = shift @ix;
            my @queue;
            for my $current (@current) {
                # print STDERR "ix: $ix, current: $current\n";
                my $dx = $ax->[$ix] - $ax->[$current];
                my $dy = $ay->[$ix] - $ay->[$current];
                if ($dx * $dx + $dy * $dy <= $radius2) {
                    # print STDERR "they are together\n";
                    push @queue, $ix;
                    last;
                }
            }
            if (@queue) {
                while (defined($ix = shift @queue)) {
                    for my $done (@done) {
                        next unless defined $done;
                        # print STDERR "ix: $ix, done: $done\n";
                        my $dx = $ax->[$ix] - $ax->[$done];
                        my $dy = $ay->[$ix] - $ay->[$done];
                        if ($dx * $dx + $dy * $dy <= $radius2) {
                            # print STDERR "they are together\n";
                            push @queue, $done;
                            undef $done;
                        }
                    }
                    push @current, $ix;
                }
            }
            else {
                push @done, $ix;
            }
        }
        if (@current >= $self->{minimum_size}) {
            @current = sort { $a <=> $b } @current if $self->{ordered};
            push @clusters, \@current;
        }
        @ix = grep defined($_), @done;
    }
    @clusters = sort { $a->[0] <=> $b->[0] } @clusters if $self->{ordered};
    return @clusters;
}

sub brute_force_clusters_ix_any {
    @_ == 1 or croak 'Usage: $clp->brute_force_clusters_ix_any';
    my $self = shift;

    my $radius = $self->{radius};
    my $radius2 = $radius * $radius;

    my $dimension = $self->{dimension};
    my $coords = $self->{coords};
    my @ix = 0..$#{$coords->[0]};
    @ix or croak "points have not been added";

    my @clusters;
    while (@ix) {
        # print STDERR "\@ix 1: ".join("-", @ix).".\n";
        my @current = shift @ix;
        my @done;
        while (@ix) {
            # print STDERR "\@ix 2: ".join("-", @ix).".\n";
            my $ix = shift @ix;
            my @queue;
            for my $current (@current) {
                # print STDERR "ix: $ix, current: $current\n";
                my $sum = 0;
                for (@$coords) {
                    my $delta = $_->[$ix] - $_->[$current];
                    $sum += $delta * $delta;
                }
                if ($sum <= $radius2) {
                    # print STDERR "they are together\n";
                    push @queue, $ix;
                    last;
                }
            }
            if (@queue) {
                while (defined($ix = shift @queue)) {
                    for my $done (@done) {
                        next unless defined $done;
                        # print STDERR "ix: $ix, done: $done\n";
                        my $sum = 0;
                        for (@$coords) {
                            my $delta = $_->[$ix] - $_->[$done];
                            $sum += $delta * $delta;
                        }
                        if ($sum <= $radius2) {
                            # print STDERR "they are together\n";
                            push @queue, $done;
                            undef $done;
                        }
                    }
                    push @current, $ix;
                }
            }
            else {
                push @done, $ix;
            }
        }
        if (@current >= $self->{minimum_size}) {
            @current = sort { $a <=> $b } @current if $self->{ordered};
            push @clusters, \@current;
        }
        @ix = grep defined($_), @done;
    }
    @clusters = sort { $a->[0] <=> $b->[0] } @clusters if $self->{ordered};
    return @clusters;
}

sub distance {
    @_ == 3 or croak 'Usage: $clp->distance($ix0, $ix1)';
    my ($self, $ix0, $ix1) = @_;
    my $coords = $self->{coords};
    my $sum = 0;
    for my $coord (@$coords) {
        my $delta = $coord->[$ix0] - $coord->[$ix1];
        $sum += $delta * $delta;
    }
    sqrt($sum);
}
1;
__END__

=head1 NAME

Algorithm::ClusterPoints - find clusters inside a set of points

=head1 SYNOPSIS

  use Algorithm::ClusterPoints;
  my $clp = Algorithm::ClusterPoints->new( dimension => 3,
                                           radius => 1.0,
                                           minimum_size => 2,
                                           ordered => 1 );
  for my $p (@points) {
      $clp->add_point($p->{x}, $p->{y}, $p->{z});
  }
  my @clusters = $clp->clusters_ix;
  for my $i (0..$#clusters) {
      print( join( ' ',
                   "cluster $i:",
                   map {
                       my ($x, $y, $z) = $clp->point_coords($_);
                       "($_: $x, $y, $z)"
                   } @{$clusters[$i]}
                 ), "\n"
           );
  }

=head1 DESCRIPTION

This module implements and algorithm to find clusters of points
inside a set.

Points can have any dimension from one to infinitum, though the
algorithm performance degrades quickly as the dimension increases (it
has O((2*D)^D) complexity).

The algorithm input parameters are:

=over 4

=item $dimension

Dimension of the problem space. For instance, for finding clusters on a
geometric plane, dimension will be 2.

=item $radius

A point P is part of a cluster when there is at least another point
from the cluster inside the (hyper-)circunference with center P and
radius $radius.

=item $minimum_size

Minimum_number of points required to form a cluster, the default is
one.

=item @points

The coordinates of the points

=item $ordered

Order the points inside the clusters by their indexes and also order
the clusters by the index of the contained points.

Ordering the output data is optinal because it can be an computational
expensive operation.

=back

=head2 API

This module has an object oriented interface with the following
methods:

=over 4

=item Algorithm::ClusterPoints->new(%args)

returns a new object.

The accepted arguments are:

=over 4

=item dimension => $dimension

number of dimensions of the points (defaul is 2).

=item radius => $radius

maximum aceptable distance between two points to form a cluster
(default is 1.0).

=item minimum_size => $minimum_size

minimun cluster size (default is 1).

=item ordered => $ordered

sort the returned data structures (default is false).

=back

=item $clp->add_point($x, $y, $z, ...)

=item $clp->add_points($x0, $y0, $z0..., $x1, $y1, $z1..., ...);

These methods register points into the algorithm.

They return the index of the (first) point added.

=item $clp->radius

=item $clp->radius($radius)

=item $clp->minimum_size

=item $clp->minimum_size($minimum_size)

=item $clp->ordered

=item $clp->ordered($ordered)

These methods get or set the algorithm parameters.

=item @coords = $clp->point_coords($index)

returns the coordinates of the point at index C<$index>.

=item @clusters_ix = $clp->clusters_ix

returns a list of clusters defined by the indexes of the points inside

The data estructure returned is a list of arrays. Every array
represents a cluster and contains the indexes of the points inside.

For instance:

  @clusters_ix = ( [ 0, 1, 5, 10, 13, 15, 17, 31, 32, 38 ],
                   [ 2, 12, 20, 26, 27, 29, 33 ],
                   [ 3, 22, 39 ],
                   [ 4, 11, 16, 30, 36 ],
                   [ 6, 14 ],
                   [ 7, 23, 24 ],
                   [ 18, 25 ],
                   [ 21, 40 ] );

You can get back the coordinates of the points using the method
C<point_coords>, as for instance:

   for my $c (@clusters_ix) {
     for my $index (@$c) {
       my ($x, $y, $z) = $clp->point_coords($index);
       ...

Or you can use the method C<clusters> described below that already
returns 2D coordinates.

=item @clusters = $clp->clusters

returns a list of clusters defined by the 2D coordinates of the points
inside.

This method is similar to C<clusters_ix> but instead of the point
indexes, it includes the point coordinates inside the cluster arrays.

This is a sample of the returned structure:

  @clusters = ( [ 0.49, 0.32, 0.55, 0.32, 0.66, 0.33 ],
                [ 0.95, 0.20, 0.83, 0.27, 0.90, 0.20 ],
                [ 0.09, 0.09, 0.01, 0.08, 0.12, 0.15 ],
                [ 0.72, 0.42, 0.67, 0.47 ],
                [ 0.83, 0.11, 0.77, 0.13, 0.73, 0.07 ],
                [ 0.37, 0.38, 0.36, 0.44 ],
                [ 0.16, 0.79, 0.14, 0.74 ] );

Note that X and Y coordinate values are interleaved inside the arrays.

=back

=head1 SEE ALSO

All began on this PerlMonks discussion:
L<http://perlmonks.org/?node_id=694892>.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2008 by Salvador FandiE<ntilde>o (sfandino@yahoo.com)

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.


=cut
