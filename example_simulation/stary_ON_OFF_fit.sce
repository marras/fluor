// ON-OFF analysis

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


//#################################
//                                #
//     On, off times division     #
//                                #
//#################################\t % 16.8f \t %16.8f \t %16.8f \t %16.8f \n', k_on_fit/k_on_teoret, sigma_on_fit/k_on_teoret, chisq_on_fit ,k_off_fit/k_off_teoret, sigma_off_fit/k_off_teoret, chisq_off_fit ) // write a line

clear{};
make_analysis = 0;
 //read theoretical parameters from file
// par_theor = fscanfMat('data/par_fit.txt')
 //bin_size = par_theor(4,1)

for(i_d = 1:1:3)
  for (bL = 1:1:size(binLs,2))
  if (WHICH_DATA(i_d) == 1)
//reading the simulation params, BUT I_on ==I_1, I_off == I_2, not always on==the higher!!
// on == on for donor!       
        if (i_d ==1)
         // I_on = par_theor(4,i_d) 
        //  I_off = par_theor(1,i_d)
          k_on_teoret = par_theor(1,i_d)
          k_off_teoret = par_theor(1,i_d)
        end
        if (i_d ==2)
        //  I_on = par_theor(4,i_d) 
        //  I_off = par_theor(1,i_d) 
          k_on_teoret = par_theor(1,i_d)
          k_off_teoret = par_theor(1,i_d)
        end
        if (i_d ==3)
        //  I_on = par_theor(4,2)/(par_theor(4,2)+par_theor(4,1)) 
        //  I_off = par_theor(1,2)/(par_theor(1,2)+par_theor(1,1))
          k_on_teoret = par_theor(1,2)
          k_off_teoret = par_theor(1,2)
          onoffANALYSIS(4) = 0;
        end

//data_in = 'data/' + WHICH_DATA(i_d) + '_1_binL_' + string(binLs(bL)) + '.txt';
        data_in = INPUT_DATA(i_d) + '_1_binL_' + string(binLs(bL)) + '.txt';

        data = fscanfMat(data_in);
        time = data(:,1); events = data(:,2);
        binL = binLs(bL)
        k_on_ini = k_on_teoret
        k_off_ini = k_off_teoret
        //thresholds
         [cutoff, cutoff_high_3, cutoff_low_3, cutoff_high_5, cutoff_low_5] = calc_cutoffs( events );
        for(thresh = 1:1:4)  
            if (thresh == 1 & onoffANALYSIS(thresh) == 1)
                thresh_u = cutoff;
                thresh_d = cutoff;
                make_analysis = 1;
            end
  
            if (thresh == 2 & onoffANALYSIS(thresh) == 1) 
                thresh_u = cutoff_high_3;
                thresh_d = cutoff_low_3;
                make_analysis = 1;
            end

            if (thresh == 3 & onoffANALYSIS(thresh) == 1)
                thresh_u = cutoff_high_5
                thresh_d = cutoff_low_5
                make_analysis = 1;
            end
  
//            if (thresh == 4 & onoffANALYSIS(thresh) == 1)
//                if (I_on > I_off)
//                    thresh_u = I_on - I_on^(0.5)
//                    thresh_d = I_off + I_off^(0.5)
//                elseif(I_off > I_on)
//                    thresh_u = I_off - I_off^(0.5)
//                    thresh_d = I_on + I_on^(0.5)
//                end
//                make_analysis = 1;
//            end
//            
            if (make_analysis == 1)
                events_size = size(events);
                states = [1:events_size(1)];

                //initial state
                if (events(1) > cutoff) then state = 1;
                elseif (events(1) <= cutoff) then state = 2;
                end //if

                //"on" and "off" states
                for i = 1:events_size(1)//*
                    if (events(i) > thresh_u) then state = 1;
                    elseif (events(i) < thresh_d) then state = 2;
                    end //if
                    states(i) = state;
                end //for*
 
                //making states data
                threshold_now = my_thresholds(thresh);
                input_data_now =  WHICH_DATA_for_file_name(i_d)
                my_times_on = 'ON_OFF_times/ON_times_' + input_data_now + '_binL_' +  string(binLs(bL)) + '_' +  threshold_now + '.txt';
                my_times_off = 'ON_OFF_times/OFF_times_' + input_data_now + '_binL_' +  string(binLs(bL)) + '_' +  threshold_now + '.txt';
                fd_on = mopen(my_times_on,'w')
                fd_off = mopen(my_times_off,'w')
                states_size = size(states);
                volume = 1;
                k=1;
                l=1;
                on = zeros();
                off = zeros();
                for j = 1:events_size(1)-1
                    if(states(j) == states(j+1)) then volume = volume + 1;
                    else time_v = volume * binL;
                        if(states(j) == 1) then on(k) = time_v;, k = k + 1;
                        elseif(states(j) == 2) then  off(l) =  time_v;, l = l+1;
                        end
                        volume = 1;
                    end
                end
                mfprintf(fd_on,'%3.9f\n',on)
                mfprintf(fd_off,'%3.9f\n',off)
                mclose(fd_on);
                mclose(fd_off);
                

//#################################
//                                #
//              Fit               #
//                                #
//#################################

                n = 15 //number of histograms bins
                //ON
                 which = 1;
                 [on_err, on_occ, on_x_mid, on_h]=histogram(on, binL, n)
// write a line

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

                 

                //cumulative on histogram  
                //[on_cum_hist] = cumulative_histogram(on_occ)
                //cum_Z_on = [on_cum_hist';1:length(on_cum_hist)'];
                // [ cum_k_on_fit, cum_sigma_on_fit, cum_chisq_on_fit] = histfit( LSW, cum_FF, dFF, cum_Z_on, k_on_ini, edges_on, hist_cutoff ); 
             
                 
                 //OFF
                 which = 2;
                 [off_err, off_occ, off_x_mid, off_h]=histogram(off, binL, n);
                 Z_off=[off_occ';1:length(off_occ)'];
                 range_off = n;
                 edges_off = (0:range_off) * off_h;
                 hist_cutoff = 0;
                 //fit
                 [ k_off_fit, sigma_off_fit, chisq_off_fit] = histfit( LSW, FF, dFF, Z_off, k_off_ini, edges_off, hist_cutoff );
                 tau_fit = 1/k_off_fit;
                 sample_size = sum( off_occ );
                 edges_left = edges_off( 1:$-1 );
                 edges_right = edges_off( 2:$ );
                 sizes = edges_right - edges_left;

                 //cumulative off histogram  
                [off_cum_hist] = cumulative_histogram(off_occ)
                 
                  
                 //writing to the file
                 my_fit = 'fit_results/ON-OFF_fit_' + input_data_now + '_binL_' + string(binLs(bL)) + '_'+ threshold_now + '.txt';
                 fd_ON_OFF = mopen(my_fit, 'a')
                 mfprintf(fd_ON_OFF, '%16.8f \t %16.8f \t % 16.8f \t %16.8f \t %16.8f \t %16.8f \n', k_on_fit/k_on_teoret, sigma_on_fit/k_on_teoret, chisq_on_fit ,k_off_fit/k_off_teoret, sigma_off_fit/k_off_teoret, chisq_off_fit ) // write a line
                 mclose(fd_ON_OFF);
              
                 
                 
//#################################
//                                #
//           Plots                #
//                                #
//#################################
                 
                if PLOTS_FACTOR == 1
                  //zmienic trzy kolejne linie (sa 2 razy)
                    edges_left = edges_on( 1:$-1 );
                    edges_right = edges_on( 2:$ );
                    sizes = edges_right - edges_left;
                    scf(3); clf(3);
                    subplot(2,1,1)
                    xtitle ("ON times fit")
                    plot2d(on_x_mid, on_occ, style = -2)
                    hist_on_avg = FF( 1:length(on_occ), k_on_teoret );//theoretical curve
                    plot2d( on_x_mid, hist_on_avg, 3 ) 
                    hist_on_fit = FF( 1:length(on_occ), k_on_fit );//fit curve
                    plot2d( on_x_mid, hist_on_fit, 5 ) 
                    subplot(2,1,2)
                    xtitle ("ON times cumulative histogram")
                    plot2d(on_x_mid, on_cum_hist, style = -4)
                    cum_hist_on_avg = cumulative_histogram(hist_on_avg');//theoretical curve
                    plot2d( on_x_mid, cum_hist_on_avg, 3 ) 

                    scf(4); clf(4);
                  //zmienic trzy kolejne linie (sa 2 razy)
                    edges_left = edges_off( 1:$-1 );
                    edges_right = edges_off( 2:$ );
                    sizes = edges_right - edges_left;
                    subplot(2,1,1)
                    xtitle ("OFF times fit")
                    plot2d(off_x_mid, off_occ, style = -2)
                    hist_off_avg = FF( 1:length(off_occ), k_off_teoret );
                    plot2d( off_x_mid, hist_off_avg, 3 )
                    hist_off_fit = FF( 1:length(off_occ), k_off_fit );
                    plot2d( off_x_mid, hist_off_fit, 5 )
                    subplot(2,1,2)
                    plot2d(off_x_mid, off_cum_hist, style = -4)
                    cum_hist_off_avg = cumulative_histogram(hist_off_avg');//theoretical curve
                    plot2d( off_x_mid, cum_hist_off_avg, 3 )                     
                    

                    //saving plots
                    filename = 'fit_results/ON_analysis_' + input_data_now + '_' + threshold_now;
                    xs2ps(3,filename)
                    unix_s(pathconvert('SCI/bin/BEpsf',%f)+' -p ' +filename)
                    filename = 'fit_results/OFF_analysis_' + input_data_now + '_' + threshold_now;
                    xs2ps(4,filename)
                    xtitle ("OFF times cumulative histogram")
                    unix_s(pathconvert('SCI/bin/BEpsf',%f)+' -p ' +filename)
                 end
                 
                 make_analysis = 0; 
               end//if (make_analysis == 1)
            end//for(thresh = 1:1:4)
            //onoffANALYSIS(4) = 1;
          end // if (WHICH_DATA(i_d) == 1)
       end // for (bL = 1:1:size(binLs,2))
    end //for(i_d = 1:1:3)
// pvm_exit;

