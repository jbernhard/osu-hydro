#!/usr/bin/perl

use Getopt::Long;

# Set defaults
#$pdgfile = "mass_width_2008-FIXED.csv";
$pdgfile = "URQMD_particles.csv";
$pbmfile = "part_tab.dat";
$pi      = 3.14159;
$pisq    = $pi*$pi;
$twopisq = 2*$pi*$pi;
$mlo     = 0;
$mhi     = 2500;
$Tstep   = 10;
$Tmax    = 200;
$Tc      = 1000;
$Tc4     = $Tc*$Tc*$Tc*$Tc;

GetOptions ( "help"    => \$help,
	     "mlo=f"   => \$mlo,
	     "mhi=f"   => \$mhi,
	     "Tstep=f" => \$Tstep,
	     "Tmax=f"  => \$Tmax,
	     "Tc=f"    => \$Tc,
	     "pbm"     => \$pbm,
	     "vh2"     => \$vh2,
	     "gev"     => \$gev,
	     "dmp"    => \$dmp
	     );

if ($help) {
    print(
	  "hrg.pl calculates EOS for Hadron Resonance Gas according to Hagedorn's prescription\n",
	  "using Modified Bessel Fn. formulas from Numerical Recipies converted to perl.\n",
	  "Writes out Temperature, Trace Anomaly, Energy Density, and Pressure (divided by T^4).\n",
	  "Usage: $0 [options]\n",
	  "Options:\n",
	  "        --help     Print this message.\n",
	  "        --mlo=f    Low  mass resonance cutoff (0 MeV) \n",
	  "        --mhi=f    High mass resonance cutoff (2500 MeV) \n",
	  "        --Tstep=f  Temperature step size in MeV (10 MeV)\n",
	  "        --Tmax=f   Maximum temperature for calculation (200 MeV)\n",
          "        --Tc=f     Match pressure at Tc to continue with 1st order transition to Bag EOS (1000 MeV)\n",
	  "        --pbm      Use part_tab.dat for resonances (default is pdg)\n",
	  "        --vh2      Use vh2_*.dat for resonances (default is pdg)\n",
	  "        --gev      Convert output temperature to GeV (input units remain MeV)\n",
	  "        --dmp      Print out mass, spin, degeneracy (m,d,s)\n"
	  );
    exit;
}

if ($pbm && $vh2) {
  print "pbm anv vh2 are incompatible options.  Please select only one.\n";
  exit;
}

# Load resonances from part_tab.dat if pbm set and form name from mass
if ($pbm) {
  open (PBM, "<$pbmfile");
  $line = <PBM>;
#  print("First line of PBMfile = $line");
  while (<PBM>) {
    ($d,$m,$s,$dummy) = split(" ");
# Convert to MeV, and multiply by -1 pre-factor
    $m *= 1000;
    $s = -$s;
    if ( ($m>mlo) && ($m<$mhi) ) {
#      print ("$m $d $s\n");
      push @m, $m;
      push @d, $d;
      push @s, $s;
      push @n, "PBM".$m;
    }
  }
  close(PBM);
}

elsif ($vh2) {
  open (VH2mass, "<vh2_pasim.dat");
  open (VH2deg,  "<vh2_gslist.dat");
  open (VH2name, "<vh2_pasinames.dat");
  while (<VH2mass>) {
    chomp;
    $m = 1000*$_;
    $n = <VH2name>;
    chomp($n);
    $d = <VH2deg>;
    chomp($d);
    # if d is even, then ferion and s=-1, else s=+1.
    $s = (2*int($d/2)==$d) ? -1 : +1;
#    print ("$n $m $d $s\n");
    push @m, $m;
    push @d, $d;
    push @s, $s;
    push @n, $n;
  }
}

# Load resonances from pdg file
else {
  open (PDG, "<$pdgfile");
  while (<PDG>) {
    if (substr($_,0,1)ne"*") {
      ($m,$mep,$men,$w,$wep,$wen,$i,$g,$j,$p,$c,$a,$mc,$q,$r,$status,$name,$quarks) = split(",");
      # convert total spin $j to float $jf and calculate $s spin symmetry
      if ($j=~/(\d+)\/2/) {
	$jfloat = $1/2.;
	$s = -1;
      } elsif ($j=~/(\d+)/) {
	$jfloat = $1;
	$s = +1;
      } else {
	$jfloat = 0;
	$s = +1;
      }
      #    print("Spin $j = $jfloat and spin_symmetry = $s\n");

      # EW and quarks have $mc<99, K0S=310, K0L=130
#      if (  {
      if ( (($mc>0)&&($mc<99)) || ($mc==310) || ($mc==130) ) {
#	print ("$name  $m $mc $q\n");
      } elsif ( ($m>mlo) && ($m<$mhi) && (($status eq "R")||($status eq "D")) ) {
	# push $m=mass, $d=degeneracy, and $j=spin onto arrays
	$dpart = (($a eq "B")||($a eq "F")) ? 2 : 1;
	$dspin = (2*$jfloat)+1;
	$d = $dpart * $dspin;
#	print ("$name  $m $j $a $mc $q $d $s\n");
	push @m, $m;
	push @d, $d;
	push @s, $s;
	push @n, $name;
      }
    }
  }
  close(PDG);
}
#Uncomment next line to test file parsing
#exit;

# Use the following to write out pdg table in array form (for root)
if ($dmp) {
  open (DMP,">mds.txt");
  $im = $#m;
  print DMP "Arrays have $im elements\n";
  print DMP "mass array\n";
  for ($i=0; $i<$#m; $i++) {$mGeV = @m[$i]/1000;print DMP "$mGeV,";}
   print DMP "\n\ndegeneracy array\n";
  for ($i=0; $i<$#d; $i++) {print DMP "@d[$i],";}
   print DMP "\n\nspin array\n";
  for ($i=0; $i<$#s; $i++) {print DMP "@s[$i],";}
  print DMP "\n";
  close (DMP);
}


#Print out lables and output format
print ("# hrg.pl with $mlo < pdg_mass < $mhi\n");
if ($gev) {
  print ("# T(GeV)   (e-3p)/T^4  e/T^4    p/T^4     Cs2 \n");
}
else {
  print ("# T(MeV)   e/T^4    p/T^4     Cs2 \n");
  print ("# T(GeV)   (e-3p)/T^4  e/T^4    p/T^4     Cs2 \n");
}
print ("# ------   ----------  --------  -----    -----\n");

$psum = 0;
$psumOld = 0;
$epsOld  = 0;
# Outer loop over T in steps of 1 MeV
for ($T=$Tstep;$T<=$Tmax;$T+=$Tstep) {

  if ($T <= $Tc) {
    # Loop over masses to calculate trace anomaly
    $msum = 0;
    for ($im=0;$im<=$#m;$im++) {
      #    print ("Working on @n[$im] @m[$im]\n");
      $ksum = 0;
      for ($ik=1;$ik<10;$ik++) {
	# if $s=1  or $ik+1 is even then sign is positive
	if ( (@s[$im]==1) ||
	     (int(($ik+1)/2.)==($ik+1)/2) ) {
	  $sign = 1;
	} else {
	  $sign = -1;
	}
	$x = @m[$im]/$T;
	$y = $ik*$x;
	$ksum += $sign*$x*$x*$x*bessk1($y)/$ik;
	#	 print ("ksumcheck : @n[$im] $ksum $sign $x $y $ik\n");
      }
      $msum += @d[$im]*$ksum/$twopisq;
      #       print ("msumcheck : @n[$im] @d[$im] $ksum $twopisq\n");
    }
    # TraceAnomaly/T^4
    $em3p = $msum;
    # Pressure/T^4, store pressure in $pc for when we cross Tc
    $psum += $em3p * $Tstep/$T;
    $pc = $psum;
    # EnergyDensity/T^4
    $eps = $em3p + 3*$psum;

    # Cs2 = delta-p / delta-E
    $T4   = $T*$T*$T*$T;
    $epsNew = $eps*$T4;
    $psumNew = $psum*$T4;
    $cs2 = ($psumNew-$psumOld)/($epsNew-$epsOld);
    $epsOld = $epsNew;
    $psumOld = $psumNew;
  }

  # continue with BAG EOS when $T exceeds $Tc (default at 1000 MeV to keep it out of the way)
  else {
#    print (">Tc=$Tc:");
    $x = $Tc/$T;
    $x4 = $x*$x*$x*$x;
    $T4   = $T*$T*$T*$T;
    $em3p = 4*(37*$pisq/90)*$x4 - 4*$pc/$T4;
    $psum += $em3p * $Tstep/$T;
    $eps = $em3p + 3*$psum;

    # Cs2 = delta-p / delta-E
    $epsNew = $eps*$T4;
    $psumNew = $psum*$T4;
    $cs2 = ($psumNew-$psumOld)/($epsNew-$epsOld);
    $epsOld = $epsNew;
    $psumOld = $psumNew;
  }

  if ($gev) {
    $Tout = $T/1000;
  }
  else {
    $Tout = $T;
  }

#  printf (" % 6.3f   % 6.4f   % 6.4f   % 6.4f\n",$T,$em3p,$psum,$eps);
  printf (" % 6.3f   % 6.4f   % 6.4f   % 6.4f   % 6.4f\n",$Tout,$em3p,$eps,$psum,$cs2);

}

exit 0;


sub bessk1 {
  local($x) = @_[0];
  local($y);
  local($ans);
  if ($x <= 2.0) {
    $y=$x*$x/4.0;
    $ans=(log($x/2.0)*bessi1($x))+(1.0/$x)*(1.0+$y*(0.15443144
	+$y*(-0.67278579+$y*(-0.18156897+$y*(-0.1919402e-1
	+$y*(-0.110404e-2+$y*(-0.4686e-4)))))));
  }
  else {
    $y=2.0/$x;
    $ans=(exp(-$x)/sqrt($x))*(1.25331414+$y*(0.23498619
	+$y*(-0.3655620e-1+$y*(0.1504268e-1+$y*(-0.780353e-2
	+$y*(0.325614e-2+$y*(-0.68245e-3)))))));
  }
}

#
# The following are perlized versions of their eponymous numerical recipes routines in C
#
sub bessk1 {
  local($x) = @_[0];
  local($y);
  local($ans);
  if ($x <= 2.0) {
    $y=$x*$x/4.0;
    $ans=(log($x/2.0)*bessi1($x))+(1.0/$x)*(1.0+$y*(0.15443144
	+$y*(-0.67278579+$y*(-0.18156897+$y*(-0.1919402e-1
	+$y*(-0.110404e-2+$y*(-0.4686e-4)))))));
  }
  else {
    $y=2.0/$x;
    $ans=(exp(-$x)/sqrt($x))*(1.25331414+$y*(0.23498619
	+$y*(-0.3655620e-1+$y*(0.1504268e-1+$y*(-0.780353e-2
	+$y*(0.325614e-2+$y*(-0.68245e-3)))))));
  }
}

sub bessi1 {
  local($x)=@_[0];
  local($ax);
  local($y);
  local($ans);

  if (($ax=abs($x)) < 3.75) {
    $y=$x/3.75;
    $y*=$y;
    $ans=$ax*(0.5+$y*(0.87890594+$y*(0.51498869+$y*(0.15084934
	 +$y*(0.2658733e-1+$y*(0.301532e-2+$y*0.32411e-3))))));
  }
  else {
    $y=3.75/$ax;
    $ans=0.2282967e-1+$y*(-0.2895312e-1+$y*(0.1787654e-1-$y*0.420059e-2));
    $ans=0.39894228+$y*(-0.3988024e-1+$y*(-0.362018e-2
	+$y*(0.163801e-2+$y*(-0.1031555e-1+$y*$ans))));
    $ans *= (exp($ax)/sqrt($ax));
  }
  return $x < 0.0 ? -$ans : $ans;
}
