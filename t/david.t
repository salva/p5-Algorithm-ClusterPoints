#!/usr/bin/perl

use strict;
use warnings;

use Test::More tests => 1;

use Algorithm::ClusterPoints;

# my @x = (0.616027832031250, 0.152069091796875, 0.190399169921875, 0.61352539062500, 0.560211181640625);
# my @y = (0.151336669921875, 0.058563232421875, 0.240966796875000, 0.65527343750000, 0.620880126953125);
# my @z = (0.325408935546875, 0.380126953125000, 0.726531982421875, 0.72015380859375, 0.141204833984375);


my @x = ( 0.377970428889331, 0.757474348942313, 0.961495921381545,
          0.903216480551723, 0.593799354197575, 0.70919235489777,
          0.0046957917138748, 0.638839705133567, 0.633172164373651,
          0.714415965376023, 0.00873674287852921, 0.161339173272754,
          0.693152409072116, 0.862988963206195, 0.139395764531077,
          0.0907198474005533, 0.142066930379361, 0.119346794697211,
          0.777107704341194, 0.890120929264217, 0.858568553548469,
          0.841220984005478, 0.976297734118152, 0.390043955935287,
          0.120808641015334, 0.322086223283808, 0.750214304217554,
          0.517964196599142, 0.909000759696017, 0.0331313886115971,
          0.669543836229327, 0.534106600125799, 0.734656711565368,
          0.490535731706249, 0.513088403946956, 0.324167543421968,
          0.359776177516881, 0.923831497415851, 0.958787652617975,
          0.549763630564346, 0.0192712245038145, 0.171134216993099,
          0.553366206475705, 0.707388792714173, 0.0099069259615483,
          0.477021927150478, 0.44888781455224, 0.396660938813032,
          0.759039345511201, 0.11884456489668, 0.644497484772689,
          0.348639407059895, 0.736432133944493, 0.0372447881542861,
          0.721632856330661, 0.0149944671694335, 0.822432989489858,
          0.89667147215869, 0.310986060345005, 0.768810123411559,
          0.995791231359807, 0.251667131638275, 0.779497318797798,
          0.0178014684322285, 0.454397736579541, 0.744710518775836,
          0.875661531753583, 0.780727837855498, 0.532228490096816,
          0.511395646951261, 0.321429102859856, 0.527783015666817,
          0.025772998173057, 0.218728490265452, 0.588982253128808,
          0.335743740377584, 0.793228815055716, 0.412322608845503,
          0.35066807300031, 0.658657761970947, 0.739743432931892,
          0.0982749499311275, 0.340369562686153, 0.26485269962885,
          0.551984785163189, 0.735503074742329, 0.614595910795533,
          0.141046950240341, 0.795321612801832, 0.0643250611283328,
          0.451573261241794, 0.414343897691115, 0.23689494438803,
          0.447063331997757, 0.376617218430415, 0.59289067293432,
          0.50786193902859, 0.646453659530042, 0.0520683407282441,
          0.675982714842608 );


my @y = ( 0.911752965084659, 0.172564697139109, 0.223033108485296,
          0.362602711067098, 0.0685045962421924, 0.418226364574053,
          0.879597770293159, 0.307493020750314, 0.60053436869471,
          0.727227930375392, 0.22964811778705, 0.120713251359575,
          0.590488382525578, 0.675663666391102, 0.0878643456820072,
          0.569853100742023, 0.298713749118249, 0.0921086046858761,
          0.123656792836172, 0.907989027480312, 0.882322350006362,
          0.662317672367777, 0.174010218679015, 0.314839684547394,
          0.544979693380085, 0.0906264865711535, 0.00463130229580599,
          0.693451138117553, 0.21507511824872, 0.478210059595796,
          0.742055614168336, 0.607568121533689, 0.0999295962595426,
          0.828355664202348, 0.449370310571894, 0.394847574430855,
          0.713697381172906, 0.83890610906397, 0.568208800628483,
          0.741070717016044, 0.226617755993587, 0.392029674374243,
          0.1996736360163, 0.250509909504217, 0.23327312775838,
          0.401472957611333, 0.446843873335453, 0.797606546710597,
          0.287509719979951, 0.9274528794584, 0.291523159129937,
          0.695116237799198, 0.80167711274693, 0.71934409988026,
          0.727884787565856, 0.963877521187612, 0.15784288048026,
          0.596318094356292, 0.125669980765892, 0.239360364608046,
          0.128638577345797, 0.326090140269599, 0.0455956418000341,
          0.100150700324605, 0.955393320323154, 0.911930463452421,
          0.979627830332547, 0.853527245506825, 0.681211219930802,
          0.448095794887791, 0.738777210925193, 0.063286857759465,
          0.586648908607415, 0.953837587567449, 0.461072723996708,
          0.354974675808535, 0.825649416749471, 0.846143143505998,
          0.462369729438521, 0.283675295343585, 0.609388646405094,
          0.961365900911318, 0.402648872886729, 0.299521230016239,
          0.535690692624158, 0.285901131511224, 0.490461924650319,
          0.0960749075949821, 0.91166826246463, 0.806063298036378,
          0.698180841295638, 0.876740451276191, 0.167244441655171,
          0.137321321859446, 0.210706299887121, 0.324361496151425,
          0.707212058913981, 0.843048010790511, 0.58029195488426,
          0.309161290830957);

my @z = ( 0.649537369389972, 0.828828198633584, 0.766584584708166,
          0.600576411217769, 0.119760777523535, 0.290968545136433,
          0.549452257176906, 0.367264137032389, 0.493365275908157,
          0.17427266706262, 0.72037474964517, 0.540358899712867,
          0.487191162691712, 0.259589772768958, 0.101884962954866,
          0.0146285473355015, 0.703503066852644, 0.0833022442925646,
          0.0484574055423899, 0.606338705080969, 0.232060756610618,
          0.28458992720682, 0.567530417864898, 0.121128637355906,
          0.294975592368413, 0.226925377508575, 0.399329666357826,
          0.76838369194574, 0.578883874185063, 0.850831315536606,
          0.325024156923657, 0.504135576198774, 0.148968333780257,
          0.964931536921323, 0.222746017951945, 0.134521715241515,
          0.0893590210699777, 0.432326530565728, 0.641005302870578,
          0.26183494222526, 0.264444835233402, 0.959325053657619,
          0.909031062228305, 0.405194770603924, 0.84562695031104,
          0.443562913713048, 0.0535699353084986, 0.877871498443124,
          0.498045385334134, 0.0885402076510324, 0.0416735112741797,
          0.253294394760513, 0.3718099689853, 0.741608267899505,
          0.119005709844284, 0.501775217897357, 0.292266948790424,
          0.162437565532525, 0.263201431501184, 0.129978315921374,
          0.793832494841805, 0.178997820438735, 0.0218077272294543,
          0.905168559060808, 0.891924348052552, 0.938081568521067,
          0.254338737785897, 0.996268354169768, 0.737522862942225,
          0.707402435927701, 0.547975652500806, 0.613754991224805,
          0.548892694953732, 0.93629937031427, 0.326442937610491,
          0.480429571852749, 0.464996382301475, 0.532078079137317,
          0.921592084121041, 0.715088769889785, 0.21513107235074,
          0.559020924979006, 0.127224125203352, 0.961130673971791,
          0.749319642528345, 0.817096527903072, 0.00468217311938801,
          0.842125962155915, 0.223987465764534, 0.00740032484203468,
          0.959745225679708, 0.0763856585734111, 0.98745607046612,
          0.441357336015653, 0.360354492232229, 0.832556722113718,
          0.459784945020882, 0.50071838838063, 0.214703489182252,
          0.158276717139795);

my @t = ( 0.669069047666817, 0.246009321367247, 0.669759760913358,
          0.672736139603423, 0.281142079463955, 0.484441321894913,
          0.694459591806456, 0.431622015711731, 0.0083859177655583,
          0.343305435170528, 0.535902828209906, 0.955785672910885,
          0.483028378492115, 0.205144729806086, 0.494054479449822,
          0.509955467586487, 0.00424239229590029, 0.907478251953599,
          0.818463696225955, 0.774827116263033, 0.843063340071232,
          0.265995338354408, 0.438260451825066, 0.526294385606356,
          0.429691642370976, 0.207674245794756, 0.420240602240561,
          0.801557223832312, 0.677086852729296, 0.830231346967601,
          0.533945575198366, 0.993513583521725, 0.944759504890662,
          0.954307142740682, 0.552642938089058, 0.56484118584158,
          0.396666963945101, 0.984795340830836, 0.16911409480624,
          0.323175476172459, 0.757119531369703, 0.147721250884565,
          0.655155361491175, 0.448790021142237, 0.718603053976373,
          0.419075745964854, 0.643698671183664, 0.716185597066023,
          0.505424286765837, 0.712011040444231, 0.66245815737917,
          0.600994219792348, 0.347085259687894, 0.138804534304093,
          0.298008866492204, 0.213604365955611, 0.865292898307921,
          0.413996973840263, 0.119203849225947, 0.250913635634767,
          0.628510148379625, 0.715111178774812, 0.738835880848136,
          0.40761510641827, 0.377954854324052, 0.567357585774126,
          0.642450060373186, 0.938423628875263, 0.678952641614835,
          0.568732858373217, 0.663213538650574, 0.312648313433023,
          0.00155005306969969, 0.0219178133062954, 0.231748625359465,
          0.152517458660647, 0.484294189048221, 0.789380913243583,
          0.779835524909593, 0.896662170958997, 0.885085101555823,
          0.938093575934381, 0.926890951248389, 0.880164925656892,
          0.473007192254038, 0.990263549756385, 0.328128871574357,
          0.772579319513547, 0.865045864377002, 0.846492605671251,
          0.313309628956709, 0.565327264627264, 0.0615956785542195,
          0.697764795586416, 0.694983693861239, 0.549879873750378,
          0.184666262539558, 0.492671527123157, 0.375136294582454,
          0.978128287010442);


my $n = @x;

my $clp = Algorithm::ClusterPoints->new(dimension => 4, ordered => 1, radius => 0.1);
$clp->add_point($x[$_], $y[$_], $z[$_], $t[$_]) for 0..$n-1;

my @bfc = $clp->brute_force_clusters_ix;
my @c = $clp->clusters_ix;

# use Data::Dumper;
# print STDERR Data::Dumper->Dump([\@bfc, \@c], [qw($bfc $c)]);
# print STDERR "distance(4, 15) = ".$clp->distance(4, 15)."\n";

is_deeply(\@c, \@bfc, "simple 4d");
# is_deeply(\@bfc, \@sol, "simple 3d brute force");

