#!/usr/bin/perl -wT

use strict;
use CGI;
use LWP::Simple;
use Template;
use Data::Dumper;
use File::Spec;
use File::Temp qw/tempdir/;
use IO::Handle;
use IPC::Open3;
use Storable;

umask 077;

my $server = "http://kinherd.org";
my $progs_path = "/scratch/localB/webseismosizer/bin";

my $self_url = "http://serv01/kiwi/cgi-bin/seismograms.pl";

if ($ENV{REMOTE_ADDR} eq '136.172.75.134') {
    $self_url = "http://localhost:8083/kiwi/cgi-bin/seismograms.pl";
}

my $gfdb_path = "/home/cesca/KINHERD/gfdb";

$ENV{PATH} = join(':', $progs_path, '/bin', '/home/cesca/wc/kiwi/bin', '/home/cesca/sw/GMT4.2.0/bin', '/home/cesca/bin', '/usr/bin' );

$ENV{LD_LIBRARY_PATH} .= ":/home/cesca/sw/netcdf-3.6.2/lib";

my $pri_color = "black";
my $sec_color = "66/103/148";

my $q = new CGI;

my $formparams = $q->Vars;

#print STDERR $ENV{REMOTE_ADDR},"\n";

my $template_url = $server."/power/trac/wiki/Template/Seismograms";
my $template = get($template_url) or die "failed to get ".$template_url;
$template =~ s|"/power/trac|"$server/power/trac|g;
$template =~ s|"/trac-htdocs|"$server/trac-htdocs|g;

my $session;
my $generation = 0;
if (exists $formparams->{session} && $formparams->{session} =~ /(\d+)/) {
    $session = $1;
} else {
    $session = int(rand(100000));
}

if (exists $formparams->{generation} && $formparams->{generation} =~ /(\d+)/) {
    $generation = $1;
}

my @errors;

my $tmpdir = File::Spec->tmpdir();
my $maindir = $tmpdir."/webseismosizer";
checkdir($maindir);
my $sessiondir = $maindir."/session".$session;
checkdir($sessiondir);
my $tempdir = tempdir( DIR => $maindir, CLEANUP => 1 ) or die "cannot create tempdir\n";
chdir($tempdir);

$ENV{HOME} = $sessiondir;

if (exists $formparams->{getfile}) {
    my $fileparam = $formparams->{getfile};
    
    ($fileparam =~ /([a-z0-9.-]+)/) or error("malformed parameter");
    
    my $file = $1;
    my $path = $sessiondir."/results-".$generation."/".$file;
    
    (-f $path) or error("file not found");
    
    ($file =~ /\.([.a-z]+$)/) or error("unrecognized filename extension");
    my $ext = $1;
    my %filetypes = ('png' => 'image/png', 'pdf' => 'image/pdf', 'tar.gz' => 'application/x-gtar' );
    (exists $filetypes{$ext}) or error("unknown filename extension");
    
    print $q->header(-type=>$filetypes{$ext}, -expires=>'now');
    open FILE,"<",$path or error("can't open file");
    binmode FILE;
    my $data;
    while (read FILE, $data, 1024) {
        print $data;
    }
    close FILE;
    exit;
}



my $javascript = <<'EOF';
<script type="text/javascript">
    var groups = new Array();
    
    function togglesource(selectelem) {
        var opts = selectelem.getElementsByTagName("option");
        for (var i=0; i<opts.length; i++) {
            var group = opts[i].value;
            var elem = document.getElementById(group);
            elem.style.display = selectelem.value == group ? 'block' : 'none';
        }
    }
</script>
EOF


my $configtemplate = <<'EOF';
<div id="[% group %]" [% IF invisible %]style="display:none"[% END %]>
    <script type="text/javascript">
        groups.push("[% group %]");
    </script>
    [% IF description %]
    <p style="text-align:justify; padding:0.2em; border:1px solid #ddc; background:#ffd">[% description %]</p> [% END %]
    [% FOREACH elements %]
    <span style="min-width:10em; display:-moz-inline-box; display:inline-block; text-align:center;">[% symbol %]:</span>
    [% IF acceptor.type == 'float' %]
        <input name="[% group %].[% symbol %]" value="[% value %]" type="text" size="5" maxlength="30"> [% unit %]
    [% ELSIF acceptor.type == 'float-table' %]
        <br>
        <textarea name="[% group %].[% symbol %]" cols="20" rows="4">[% value %]</textarea>
    [% ELSIF acceptor.type == 'select' %]
        <select name="[% group %].[% symbol %]" size="1" onchange="[% acceptor.options.action %]">
            [% FOREACH  acceptor.args %]
                <option[% IF arg == value %] selected[% END %]>[% arg %]</option>
            [% END %]
        </select>
    [% END %]
    <br>
[% END %]
</div>
EOF

my $sourceinfos = source_info();
my @source_types = map { $_->{type} } @$sourceinfos;

my $defaultreceivers = << 'EOF';
42.350 13.400
49.780 17.540
45.490 25.950
EOF

#47.920 19.890
#35.870 14.520
#34.960 33.330
#35.280 24.890
#35.180 25.500
#49.630 22.710
#36.370 25.460
#42.620 23.240
#EOF
my @gfdbs_dirty = grep { -d $gfdb_path."/".$_ && -f $gfdb_path."/".$_."/db.index" } split(' ', `ls $gfdb_path` );
my @gfdbs;
foreach (@gfdbs_dirty) {
    if (/^([a-zA-Z0-9_-]+)$/) {
        push @gfdbs, $1;
    }
}


my $results = previous_results( $sessiondir );

# params printed in form should always correspond to the latest results
if (@$results) {
    my $lformp = $results->[-1]{formparams};
    foreach my $k (keys %$lformp) {
        if ($k =~ /^(source|receivers|bilateral|circular|gfdb)\./) {
            if (! exists $formparams->{$k}) {
                $formparams->{$k} = $lformp->{$k};
            }
        }
    }
}

my %description = ( 
'bilateral' =>
'This sourcetype can be used to make unilateral or bilateral Haskell type ruptures. '.
'It must be taken care of keeping the source completely inside the range of provided Greens functions. (And to keep it completely below the earth\'s surface)',
'circular' =>
'This sourcetype can be used to make circular ruptures. '.
'It must be taken care of keeping the source completely inside the range of provided Greens functions. (And to keep it completely below the earth\'s surface)' );

my %source_config = (
    'group' => 'source',
    'elements' => [
        ele( 'source', 'latitude', 'degrees', acc('float', -90.,90.), 40.75),
        ele( 'source', 'longitude', 'degrees', acc('float', -180.,180.), 29.86),
        ele( 'source', 'effective-dt', 's', acc('float', 0.1, 10.), 0.5),
        ele( 'source', 'sourcetype', undef, acc('select', @source_types, { action => "togglesource(this)" } ), $source_types[0]) ]
);
my @source_type_configs;
foreach my $sourceinfo (@$sourceinfos) {
    push @source_type_configs, {
        'group' => $sourceinfo->{type},
        'description' => $description{$sourceinfo->{type}},
        'invisible' => ($sourceinfo->{type} eq $source_config{elements}[3]{value}) ? 0 : 1,
        'elements' => [ 
            map { ele( $sourceinfo->{type}, $_->{name}, $_->{unit}, 
                        acc('float', $_->{soft_min}, $_->{soft_max} ), 1.*$_->{default}) } @{$sourceinfo->{params}} ]
    };
}

my %receivers_config = (
    'group' => 'receivers',
    'elements' => [ ele( 'receivers', 'latitude-longitude', 'degrees', acc('float-table', [-90.,90.], [-180.,180.] ), $defaultreceivers ) ]
);

my %gfdb_config = (
    'group' => 'gfdb',
    'elements' => [ ele('gfdb','greens-function-database', undef, acc('select', @gfdbs, 
        { action => "togglesource(this)" } ), 'marmak135' ) ]
);

my $tt = Template->new( { INTERPOLATE => 1 } );

my @groups = ( \%source_config, @source_type_configs, \%receivers_config, \%gfdb_config );

foreach my $group (@groups) {
    my $output = '';
    $tt->process(\$configtemplate,$group, \$output) || die $tt->error(),"\n";
    $group->{htmldata} = $output;
}


if ($formparams->{calculate}) {
    $generation = scalar @$results +1;
    my $result = calculate( $formparams ) if $formparams->{calculate};
    if ($result) {
        if ($generation > 1) {
            $result->{changes} = results_dif( $results->[-1]{formparams}, $result->{formparams} );
        }
        push @$results, $result;
    }
}

my @images;
my $plotgeneration = 'not available';
my $plotgenerationref = 'not available';
if (exists $formparams->{display}) {
    if ($formparams->{display} =~ /(seismograms|spectra)/) {
        my $displaywhat = $1;
        my $plots = plot( $generation, $generation-1, $displaywhat );
        push @images, @$plots;
        $plotgeneration = $generation;
        $plotgenerationref = $generation-1 unless $generation-1 == 0;
    }
}

my %vars = map { $_->{group} => $_->{htmldata} } @groups;

my @gfdbinfos;
foreach (@gfdbs) {
    my $info = gfdb_info( $gfdb_path."/".$_."/db",
                          $gfdb_path."/".$_."/gfdb_info_output.txt",
                          $gfdb_path."/".$_."/info.txt");
    $info->{name} = $_;
    $info->{invisible} = ($info->{name} eq $gfdb_config{elements}[0]{value}) ? 0 : 1,
    push @gfdbinfos, $info;
}

$vars{javascript} = $javascript;
$vars{errors} = join("\n", @errors);
$vars{session} = $session;
$vars{generation} = $generation;
$vars{results} = $results;
$vars{images} = \@images;
$vars{selfurl} = $self_url;
$vars{gfdbinfos} = \@gfdbinfos;
$vars{plotgeneration} = $plotgeneration;
$vars{plotgenerationref} = $plotgenerationref;

print $q->header;
$tt->process(\$template, \%vars) || die $tt->error(),"\n";

sub gfdb_info {
    my ($dbpath,$gfdb_info_output, $infotextfile) = @_;
    my %i;
    if (-f $gfdb_info_output) {
        %i = map { chomp; split('=',$_,2) } `cat $gfdb_info_output`;
    } else {
        %i = map { chomp; split('=',$_,2) } `gfdb_info $dbpath`;
    }
    my $minxkm = $i{firstx}/1000.;
    my $maxxkm = ($i{firstx}+($i{nx}-1)*$i{dx})/1000.;
    my $minzkm = $i{firstz}/1000.;
    my $maxzkm = ($i{firstz}+($i{nz}-1)*$i{dz})/1000.;
    $i{xrangekm} = $minxkm ." - ". $maxxkm;
    $i{dxkm} = $i{dx}/1000.;
    $i{zrangekm} = $minzkm ." - ". $maxzkm;
    $i{dzkm} = $i{dz}/1000.;
    
    if (-f $infotextfile) {
        $i{infotext} = `cat $infotextfile`;
    } else {
        $i{infotext} = '';
    }
    
    return \%i;
}


sub previous_results {
    my @results;
    my $igen = 1;
    while (1) {
        my $path = $sessiondir . "/results-".$igen;
        last unless (-d $path);
        my %entry;
        $entry{generation} = $igen;
        $entry{seismogram_plots} = 1 if (-f $path.'/seismogram-1-1.png');
        $entry{spectrum_plots} = 1   if (-f $path.'/spectra-1-1.png');
        $entry{formparams} =  retrieve($path.'/formparams');
        if ($igen > 1) {
            $entry{changes} = results_dif( $results[-1]{formparams}, $entry{formparams} );
        }
        push @results, \%entry;
        $igen ++;
    }
    closedir(DIR);
    return \@results;
}

sub results_dif {
    my ($a,$b) = @_;
    my @changes;
    foreach (keys %$a) {
        next if ($_ eq '__FORM_TOKEN');
        next if ($_ eq 'generation');
        if (exists $a->{$_} && exists $b->{$_} && $a->{$_} ne $b->{$_}) {
            my %change;
            $change{key} = $_;
            $change{from} = $a->{$_};
            $change{to} = $b->{$_};
            push @changes, \%change;
        }
    }
    return \@changes;
}

sub calculate {
    my ($formparams) = @_;

    my $receiverfile = $tempdir."/receivers";
    my $recs = $formparams->{'receivers.latitude-longitude'};
    chomp($recs);
    my @recs = split("\n", $recs);
    open OUT, ">".$receiverfile;
    foreach (@recs) {
        print OUT $_," ard\n";
    }
    close OUT;

    my @receivers = split("\n",$formparams->{'receivers.latitude-longitude'});
    my $nrec = scalar(@receivers);
    
    my $gfdb = $gfdb_path."/".$formparams->{'gfdb.greens-function-database'}."/db";
    my $effective_dt = $formparams->{'source.effective-dt'};
    my $sourcelat = $formparams->{'source.latitude'};
    my $sourcelon = $formparams->{'source.longitude'};
    my $sourcetype = $formparams->{'source.sourcetype'};
    my ($sourceinfo) = grep { $_->{type} eq $sourcetype } @$sourceinfos;
    my $sourceparams = join(" ", map { $formparams->{$sourcetype.".".$_->{name}} } @{$sourceinfo->{params}} );
    my $seismogramstem = 'seismogram';
    my $seismogrambase = $tempdir."/".$seismogramstem;
    my $spectrastem = 'spectrum';
    my $spectrabase = $tempdir."/".$spectrastem;
    
    
    my $commands = <<"EOF";
set_database $gfdb
set_effective_dt $effective_dt
set_receivers $receiverfile
set_source_location $sourcelat $sourcelon 0
set_source_params $sourcetype $sourceparams
calculate_seismograms
output_seismograms $seismogrambase table
output_seismogram_spectra synthetics $spectrabase
EOF

    open ERR, ">".$tempdir."/errors";
    my $null = "";
    my $pid = open3( \*TO, \*FROM, ">&ERR", 'minimizer' );

    mm($commands);
    
    close TO;
    close FROM;
    close ERR;
    
    wait;
    my $failed = 0;
    if ($? != 0) {
        push @errors, "minimizer had a bad exit status.";
        $failed = 1;
    }
    
    open ERR, "<".$tempdir."/errors";
    my $count = 0;
    while (<ERR>) {
        print STDERR $_;
        $count++;
    }
    close ERR;
    if ($count != 0) {
        push @errors, "minimizer threw up some error messages.";
        return undef;
    }
    return undef if ($failed);
    
    # package files
    
    my $packagedir = $tempdir."/results-".$generation;
    mkdir $packagedir;
    $formparams->{nrec} = $nrec;
    for (my $irec=1; $irec<=$nrec; $irec++) {
        foreach my $icomp (qw/a r d/) {
            foreach my $base ($seismogrambase,$spectrabase) {
                my $file = $base."-".$irec."-".$icomp.".table";
                ($file =~ m|([^/]+$)|) or next;
                my $filenopath = $1;
                rename $file, $packagedir."/".$filenopath;
            }
        }
    }
    # package results
    rename $receiverfile, $packagedir."/receivers.table";
    my $tarfile = "results-".$generation.".tar.gz";
    system('tar', '-czf', $tempdir."/".$tarfile, "results-".$generation );
    
    # put stuff to results directory
    my $resultsdir = $sessiondir."/results-".$generation;
    rename $packagedir, $resultsdir or die "cannot move packagedir\n";
    rename $tempdir."/".$tarfile, $resultsdir."/".$tarfile or die "cannot move tarfile\n";
    
    store $formparams, $resultsdir."/formparams";
    
    return { 'generation' => $generation, 'formparams' => $formparams };
}

sub plot {
    my ($generation, $ref_generation, $kind) = @_;
    my $resultsdir = $sessiondir."/results-".$generation;
    opendir(DIR, $resultsdir) or die "cannot opendir: ".$resultsdir."\n";
    
    my $formparams = $results->[$generation-1]->{formparams};
    my $nrec = $formparams->{nrec};
    my @receivers = split("\n",$formparams->{'receivers.latitude-longitude'});
    
    if ($ref_generation) {
        my $ref_formparams = $results->[$ref_generation-1]->{formparams};
        if ($formparams->{'receivers.latitude-longitude'} ne
            $ref_formparams->{'receivers.latitude-longitude'}) {
            $ref_generation = 0;
        }
    }
    my $ref_resultsdir = $sessiondir."/results-".$ref_generation;
    
    my %km = ( 'seismograms' => 'seismogram' , 'spectra' => 'spectrum' );
    my %plotfunc = ( 'seismograms' => \&plot_seismograms, 'spectra' => \&plot_spectra );
    my @output;
    for (my $irec=1; $irec<=$nrec; $irec++) {
        my @infiles;
        my @refiles;
        foreach my $icomp (qw/a r d/) {
            my $id = "-".($irec)."-".($icomp);
            my $infile = $resultsdir.'/'.$km{$kind}.$id.".table";
            my $refile = $ref_resultsdir.'/'.$km{$kind}.$id.".table" if ($ref_generation);
            push @infiles, $infile;
            push @refiles, $refile if $refile;
        }
        my $outbase = $resultsdir.'/'.$km{$kind}.'-'.$irec;
        my $stem = $km{$kind}.'-'.$irec;
        
        &{$plotfunc{$kind}}( \@infiles, \@refiles, $outbase);
        
        my ($wid,$hei) = pngsize( $outbase.".png" );
        
        push @output, { pngimage => $stem.".png", 
                        pdfimage => $stem.".pdf", 
                        ireceiver => $irec,
                        generation => $generation,
                        width => $wid,
                        height => $hei,
                        receiver => $receivers[$irec-1] };
    }
    return \@output;
}

sub plot_seismograms {
    my ($infiles, $refiles, $outstem ) = @_;
    
    my $pdffile = $outstem.".pdf";
    my $pngfile = $outstem.".png"; 
    return if (-f $pngfile);
    
    my ($xmin,$xmax,$ymin,$ymax) = minmaxfiles( @$infiles,@$refiles);
    my @options = qw/--ydistrib --axes=S --height=2 --width=9 --xlabel=Time --xunit=s --fit --topmargin=0.2 --xautoscale=min-max/;
    push @options, '--xrange='.$xmin."/".$xmax;
    push @options, '--yrange='.$ymin."/".$ymax;
    
    system('gmtset', 'LABEL_FONT_SIZE', '=', 14);
    if (@$refiles) {
        system('autoplot',(map {$_.'[1p/'.$sec_color.']'} reverse @$refiles), $pdffile, @options, "-K" );
        push @options, "-O";
    }
    system('autoplot',(map {$_.'[1p/'.$pri_color.']'} reverse @$infiles), $pdffile, @options);
    system('convert',$pdffile,$pngfile );
}

sub plot_spectra {
    my ($infiles, $refiles, $outstem ) = @_;
    
    my $pdffile = $outstem.".pdf";
    my $pngfile = $outstem.".png"; 
    return if (-f $pngfile);
    my ($xmin,$xmax,$ymin,$ymax) = minmaxfiles( @$infiles,@$refiles);
    
    my @options = qw/--axes=S --height=2 --width=9 --xlabel=Time 
                    --xunit=s --fit --topmargin=0.2
                    --xautoscale=min-max --xlog
                    --yautoscale=min-max --ylog/;
    push @options, '--xrange=0.01/1';
    $ymin = $ymax/10000.;
    push @options, '--yrange='.$ymin."/".$ymax;
    system('gmtset', 'LABEL_FONT_SIZE', '=', 14);
    
    my @files;
    my $ntr;
    if (@$refiles) {
        $ntr = 2;
        push @files, @$refiles;
        push @files, @$infiles;
    } else {
        $ntr = 1;
        @files = @$infiles;
    }
    multiplot( $ntr, @files, $pdffile, "--", @options );
    system('convert',$pdffile,$pngfile );
}


sub multiplot {

    my $ntr = shift @_;
    my $n = 0;
    foreach (@_) {
        last if ($_ eq '--');
        $n++;
    }
    
    $n--;
    my $outfile = $_[$n];
    
    $n % $ntr == 0 or die "$n % $ntr != 0";
    $n = $n/$ntr;
    my @files;
    for (my $itr=0; $itr<$ntr; $itr++) {
        my @filegroup = @_[$n*$itr .. $n*($itr+1)-1];
        push @files, \@filegroup;
    }
    my @addargs = @_[$n*$ntr+2 .. @_-1];
    
    my $nx = 3;
    my $ny = int(($n+1)/$nx);
    
    foreach my $i (0 .. $n-1) {
        my $ix = $i % $nx + 1;
        my $iy = int($i / $nx) + 1;
        my @ko;
        push @ko, '-K' if ($i != $n-1);
        push @ko, '-O' if ($i != 0);
        
        my @fs = map { $_->[$i] } @files;
        my @linestyles;
        if ($ntr == 1) {
            @linestyles = ("[1p/".$pri_color."]");
        } elsif ($ntr == 2) {
            @linestyles = ("[1p/".$sec_color."]","[1p/".$pri_color."]");
        } else {
            @linestyles = ("[1p/red]","[1p/green]","[1p/blue]");
        }
        foreach (0..$ntr-1) {
            $fs[$_] .= $linestyles[$_];
        }
      
        system('autoplot', @fs, $outfile, @ko,
            "--xnlayout=".$nx,
            "--ynlayout=".$ny,
            "--xilayout=".$ix,
            "--yilayout=".$iy, @addargs );
    }
}

sub pngsize {
    my ($file) = @_;
    my $fout = `file $file`;
    if ($fout =~ /(\d+) x (\d+)/) {
        return ($1,$2);
    } else {
        return (0,0);
    }
}

sub minmaxfiles {
    my @files = @_;
    my ($xmin,$xmax,$ymin,$ymax);
    foreach my $infile (@files) {
        open IN,"<",$infile or next;
        while (<IN>) {
            my ($x,$y) = split;
            $xmin = $x if (!defined($xmin) || ($x < $xmin));
            $xmax = $x if (!defined($xmax) || ($x > $xmax));
            $ymin = $y if (!defined($ymin) || ($y < $ymin));
            $ymax = $y if (!defined($ymax) || ($y > $ymax));
        }
        close IN;
    }
    my ($xcmin, $xcmax) =  map { m/([0-9.eE+-]+)/ ? $1 : 0 } ($xmin,$xmax); # launder for taint check
    my ($ycmin, $ycmax) =  map { m/([0-9.eE+-]+)/ ? $1 : 0 } ($ymin,$ymax);
    
    return ($xcmin, $xcmax, $ycmin, $ycmax);
}

# communicate with minmizer
# push command to minimizer, get answer
sub mm {
    my ($cmds) = @_;
    chomp($cmds);
    my @cmds = split(/\n/,$cmds);
    my @answers;
    foreach my $cmd (@cmds) {
        print TO $cmd,"\n";

        my $status = <FROM>;
        return undef unless defined $status;
        if ($status =~ /ok >$/) {
            my $answer = <FROM>;
            chomp($answer);
            push @answers, $answer;
        }
    }
    return @answers;
}


sub checkdir {
    my ($dir) = @_;
    (-d $dir) or mkdir $dir or die "cannot create directory: ".$dir."\n";
    (-o $dir) or die "refusing to use foreign directory: ".$dir."\n";
    my @stats = stat($dir);
    ($stats[2] & 0777) == 0700 or die "working directory has wrong access mode: ".$dir."\n";
}

sub ele {
    my ($group, $symbol, $unit, $acceptor, $default) = @_;
    my $value = $default;
    if (exists $formparams->{$group.'.'.$symbol}) {
        $value = $formparams->{$group.'.'.$symbol};
    }
    my %config = (
        group => $group,
        symbol => $symbol,
        unit => $unit,
        acceptor => $acceptor,
        default => $default,
        value => $value
    );
    my $correctedvalue = correct(\%config);
    $formparams->{$group.'.'.$symbol} = $correctedvalue;
    $config{value} = $correctedvalue;
    return \%config;
}

sub acc {
    my %acc;
    $acc{type} = $_[0];
    my @args;
    if ( scalar(@_) > 1 ) {
        foreach (@_[1..(scalar(@_)-1)]) {
            if (ref($_) eq "HASH") {
                $acc{options} = $_;
            } else {
                push @args, $_;
            }
        }
    }
    my @xx = map { {'arg' => $_} } @args;
    $acc{args} = \@xx;
    return \%acc;
}

sub correct {
    my ($ele) = @_;
    my $acc = $ele->{acceptor};
    my $value = $ele->{value};
    if ($acc->{type} eq 'float') {
        if ($value =~ /([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/) {
            my $float = $1;
            my $min = $acc->{args}[0]{arg};
            my $max = $acc->{args}[1]{arg};
            if ($float < $min || $max < $float) {
                push @errors, "Value for ".$ele->{symbol}." must be in range [".(1.*$min).",".(1.*$max)."]. Resetting it to its default value of ".$ele->{default}.".";
                return $ele->{default};
            }
            return $float;
            
        } else {
            push @errors, "Expected float value for ".$ele->{symbol}.". Resetting it to its default value of ".$ele->{default}.".";
            return $ele->{default};
       }
    } elsif ($acc->{type} eq 'select') {
        my %allowed = map { $_->{arg} => 1 } @{$acc->{args}};
        if (! exists $allowed{$value}) {
            push @errors, "Expected one the following values for ".$ele->{symbol}.": {".
                 join(", ", keys(%allowed))."}. Resetting it to its default value of ".$ele->{default}.".";
        } else {
            return $value
        }
    } elsif ($acc->{type} eq 'float-table') {
        my @lines = split(/[\n\r]+/,$value);
        my @args = map { $_->{arg} } @{$ele->{acceptor}->{args}};
        my $wantcols = scalar(@args);
        my @mins; my @maxs;
        foreach (@args) {
            push @mins, $_->[0];
            push @maxs, $_->[1];
        } 
        my $iline = 0;
        my @clean;
        LINE: foreach my $line (@lines) {
            $iline++;
            my @tokens = split(" ",$line);
            if (scalar(@tokens) != $wantcols) {
                push @errors, "Expected ".$wantcols." column".(($wantcols>1)?'s':'').
                     " in table ".$ele->{symbol}." at line ".$iline.". Deleting this line.";
                next LINE;
            }
            my $icol = 0;
            my @cleancols;
            foreach my $token (@tokens) {
                $icol++;
                my $col = '';
                if ($token =~ /([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)/) {
                    my $float = $1;
                    if ($float < $mins[$icol-1] || $maxs[$icol-1] < $float) {
                        push @errors, "Value in column ".$icol." in table ".$ele->{symbol}." at line ".$iline." must be in range [".(1.*$mins[$icol-1]).",".(1.*$maxs[$icol-1])."]. Deleting this line.";
                        next LINE;
                    } else {
                        push @cleancols, 1.*$float;
                    }
                } else {
                    push @errors, "Expected float value in column ".$icol." in table ".$ele->{symbol}." at line ".$iline.". Deleting this line.";
                    next LINE;
                }
                
            }
            push @clean, join(" ",@cleancols);
        }
        if (scalar(@clean)<1) {
            push @errors, "Need at least one row of valid data in table ".$ele->{symbol}.". Resetting table to its default.";
            return $ele->{default};
        }
        return join("\n",@clean)."\n";
    }
    push @errors, "Unable to validate value for ".$ele->{symbol}.", because no validation is implemented for it. Please inform the developer. Resetting it to its default value of '".$ele->{default}."'.";
    return $ele->{default};
}

# suck source_info output into one data structure
sub source_info {
    my $sourceinfo = `source_info`;
    $sourceinfo =~ s/^source types://;
    my @sourcetypes = split(" ", $sourceinfo);
    my @sourceinfos;
    foreach my $sourcetype (@sourcetypes) {
        if ($sourcetype =~ /([a-z]+)/) {
            my %params;
            my $sourcetype_copy = $1;
            foreach my $line (`source_info $sourcetype_copy`) {
                foreach ("parameter names:", "parameter units:", "parameter hard min:",
                        "parameter hard max:", "parameter soft min:",
                        "parameter soft max:", "parameter defaults:") {
                    if ($line =~ s/$_//) {
                        my $k = $_;
                        $k =~ s/parameter //;
                        $k =~ tr/://d;
                        $k =~ tr/ /_/;
                        $k =~ s/s$//; # singular
                        $params{$k} = [ split(' ',$line) ];
                    }
                }
            }
            my @params;
            for (my $i=0; $i<scalar(@{$params{name}}); $i++) {
                my %param;
                foreach my $k (keys %params) {
                    $param{$k} = $params{$k}->[$i];
                }
                push @params, \%param;
            }
            my %sourceinfo;
            $sourceinfo{type} = $sourcetype;
            $sourceinfo{params} = \@params;
            push @sourceinfos, \%sourceinfo;
        }
    }
    return \@sourceinfos;
}

sub error {
    print $q->header("text/plain");
    print $_[0];
    exit;
}
