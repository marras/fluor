//#################################
//                                #
//            FUNCTIONS           #
//                                #
//#################################

//fit_on-function
function [e]=G(p,z)
   if z(1) > 0 then
     e = (z(2)-size_on_off(1)*h*exp(-z(1)/p(1))/p(1))./sqrt(z(1)) 
   else
     e = 0
   end //if
endfunction

// LS fit function 
function y = FF(hist_i, p)
    y=sample_size*exp(-p(1).*(edges_left(hist_i))).*(1.0-exp(-p(1).*sizes(hist_i)))
endfunction

// derivative
function y = dFF(hist_i, p)
    y=sample_size*exp(-p(1).*(edges_left(hist_i))).*(-(edges_left(hist_i))+(edges_right(hist_i)).*exp(-p(1)*sizes(hist_i)))
endfunction

// LS weigthed cost function 
function e = LSW(p, z),
    y = z(1), x = z(2);
    if ( y > hist_cutoff ) then 
      e = (y - FF(x, p)) / sqrt(y)
      else e = 0.0;
    end
endfunction

// histogram fit
function [ k_fit, sigma_fit, chisq_fit ] = histfit(LSW, FF, dFF, Z, k_ini, edges, hist_cutoff)
    edges_left = edges( 1:$-1 );
    edges_right = edges( 2:$ );
    sizes = edges_right - edges_left;
    sample_size = sum( Z(1,:) );
    [ k_fit, err_fit ] = datafit( LSW, Z, k_ini );
    chisq_fit = err_fit / ( length( find( Z(1,:) > hist_cutoff ) ) -1 );
    
    F = zeros( length( Z(1,:) ), 1);
    for i = 1:length( Z(1,:) ), 
        if ( Z(1,i) > hist_cutoff ) F(i,1) = dFF(i, k_fit) ./ sqrt( Z(1,i) );
        else F(i,1) = 0.0;
        end
    end  

    cov_fit = inv( F'*F )
    sigma_fit = sqrt( cov_fit )
endfunction

//thresholds
function [cutoff, cutoff_high_3, cutoff_low_3, cutoff_high_5, cutoff_low_5] = calc_cutoffs( events )
   // initial estimate
   cutoff = ( max(events) + min(events) ) / 2.0; 
   while(1)
      avg_off = sum( events(find(events<cutoff)) ) / length(events(find(events<cutoff)));
      avg_on = sum( events(find(events>cutoff)) ) / length(events(find(events>cutoff)));
      cutoff_new = ( avg_off + avg_on ) / 2.0
      if abs( cutoff_new - cutoff ) / cutoff < 0.00001 then break;
      else cutoff = cutoff_new;
      end
    end;//while
    cutoff_low_3 = avg_off + (1/3.0)*(avg_on - avg_off);
    cutoff_high_3 = avg_off + (2/3.0)*(avg_on - avg_off);
    cutoff_low_5 = avg_off + (1/5.0)*(avg_on - avg_off);
    cutoff_high_5 = avg_off + (4/5.0)*(avg_on - avg_off);
endfunction

// derivative $$$ need to standard dev
function y = dFF(hist_i,p)
    y=sample_size*exp(-p(1).*(edges_left(hist_i))).*(-(edges_left(hist_i))+(edges_right(hist_i)).*exp(-p(1)*sizes(hist_i)))
endfunction

//vizualization, theoretical curve and fit
function  [err, occ, x_mid, h, n]=histogram(on_off, binL, n)
 //"x" for histograms
 maxim = max(on_off);
 minim = min(on_off);
 size_on_off = size(on_off);
 //bin size for histogram
 //n = 15; //start numbers of bins
 h_temp = (maxim - minim)/n
 h = (ceil(h_temp/binL)+1)*binL;
 x_mid=[minim+h/2:h:n*h+h/2]';
 val = linspace(minim,n*h,n+1)'
 [err, occ] = dsearch(on_off, val)
endfunction

//cumulative histograms
function [cum_hist]=cumulative_histogram(on_off_occ)
  tmp_occ = 0;
  occ_size = size(on_off_occ)
  for(ch = 1:1:occ_size(1))
    occ_sum = sum(on_off_occ)
    tmp_occ = tmp_occ + on_off_occ(ch)
    cum_hist(ch) = tmp_occ/occ_sum
  end  
endfunction
 on = fscanfMat('on_final.dat')
binL = 0.5;
k_on_ini = 1;

                n = 100 //number of histograms bins
                //ON
                 which = 1;
                 [on_err, on_occ, on_x_mid, on_h]=histogram(on, binL, n)
                 fd_h = mopen('hist_on.txt', 'w');
                 mfprintf(fd_h, '%8.2f \t %8.2f \n', on_x_mid(:,:), on_occ(:,:));
                 mclose(fd_h);
                 Z_on=[on_occ';1:length(on_occ)'];
                 range_on = n;
                 edges_on = (0:range_on) * on_h;
                 hist_cuton = 0;
                 //fit
                 [ k_on_fit, sigma_on_fit, chisq_on_fit] = histfit( LSW, FF, dFF, Z_on, k_on_ini, edges_on, hist_cuton );
                 tau_fit = 1/k_on_fit;
                 sample_size = sum( on_occ );
                 edges_left = edges_on( 1:$-1 );
                 edges_right = edges_on( 2:$ );
                 sizes = edges_right - edges_left;

scf(0); clf(0);
subplot(2,1,1);
plot2d(on_x_mid, on_occ, -2);
//hist_on_avg = FF( 1:length(on_occ), k_on_teoret );//theoretical curve
//plot2d( on_x_mid, hist_on_avg, 3 ) 
hist_on_fit = FF( 1:length(on_occ), k_on_fit );//fit curve
plot2d( on_x_mid, hist_on_fit, 5 ) 
ti = 'k_on = ' + string(k_on_fit);
xtitle(ti);

nevents = sum(on_occ);
subplot(2,1,2);
xtitle ("Differences between simulated and fitted values")
plot2d(on_x_mid, (on_occ - hist_on_fit')./sqrt(hist_on_fit'/nevents), style=-2)//./sqrt(fitHist'/nevents)

filename='ON_hist';
xs2pdf(0,filename);
xs2pdf(gcf(),filename);
