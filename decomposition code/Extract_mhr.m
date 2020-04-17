function HR = Extract_hr(rtfr_post, basicTF, dtw)


rtfr_post(1:round(0.4/basicTF.fr),:)=0;
rtfr_post(round(2.4/basicTF.fr):end,:)=0;
HR = CurveExt_M(rtfr_post(1:round(2.4/basicTF.fr),:)', dtw);

% for ci = 1:length(HR_ma)
%     low_bnd = round(HR_ma(ci)*(1-0.05)); high_bnd = round(HR_ma(ci)*(1+0.05));
%     rtfr_post(low_bnd:high_bnd,ci) = rtfr_post(low_bnd:high_bnd,ci)./10000;
%     rtfr_post(2*low_bnd:2*high_bnd,ci) = rtfr_post(2*low_bnd:2*high_bnd,ci)./10000;
%     rtfr_post(round(0.5*low_bnd):round(0.5*high_bnd),ci) = rtfr_post(round(0.5*low_bnd):round(0.5*high_bnd),ci)./10000;
% end
% 
%[HR_fe] = CurveExt_M(rtfr_post(round(0.9/basicTF.fr)+1:round(3.2/basicTF.fr),:)', dtw);
%HR_fe = HR_fe + round(0.9/basicTF.fr) ;
%
%if mean(HR_ma) > mean(HR_fe)
%	HR = HR_fe ;
%else
%    HR = HR_ma;
%end
%
%% p_HR_ma = basicTF.fs./HR_ma./basicTF.fr;
%% p_HR_fe = basicTF.fs./HR_fe./basicTF.fr;
%% 
%% z = t; %basicTF.fs*time_stamp.*(0:size(rtfr_post,2));
%% start_idx = floor(param.seg_len/2/basicTF.hop)+1;
%% for ti = start_idx:size(rtfr_post,2)-start_idx
%% %     z = basicTF.fs*time_stamp*(ti-1)+1;
%%     p = [round(p_HR_ma(ti)*0.94):round(p_HR_ma(ti)*1.06) round(p_HR_fe(ti)*0.94):round(p_HR_fe(ti)*1.06)];
%%     Ax = sum(round(p_HR_ma(ti)*0.94):round(p_HR_ma(ti)*1.06)); 
%%     Bx = sum(round(p_HR_fe(ti)*0.94):round(p_HR_fe(ti)*1.06));
%%     W = [];
%%     for pi = 1:length(p)
%%         W = [W; p(pi).^param.power_weight.*ones(p(pi),1)];
%%     end
%% 
%%     [A, Phi] = impulse_train_dict(param.seg_len, p);
%%     A = A./repmat(sqrt(sum(A.^2)),[size(A,1) 1]);
%% %     max([1 z(ti)-floor(param.seg_len/2)]:min([z(ti)+floor(param.seg_len/2) length(x)];
%%     xx = x(z(ti)-floor(param.seg_len/2):z(ti)+floor(param.seg_len/2))';
%% %     y = mexLasso(xx, A, param);%A\xx';
%%     y = mexLassoWeighted(xx, A, W, param);
%%     MM=A(:,1:Ax)*y(1:Ax);
%%     FF=A(:,Ax+1:end)*y(Ax+1:end);
%% %     mix_sig(1+(ti-200)*20:200+(ti-200)*20) = mix_sig(1+(ti-200)*20:200+(ti-200)*20)+ xx.*hann(200);
%%     ma_sig(z(ti)-floor(param.seg_len/2):z(ti)+floor(param.seg_len/2)) = ma_sig(z(ti)-floor(param.seg_len/2):z(ti)+floor(param.seg_len/2))+ (MM.*hann(param.seg_len))';
%%     fe_sig(z(ti)-floor(param.seg_len/2):z(ti)+floor(param.seg_len/2)) = fe_sig(z(ti)-floor(param.seg_len/2):z(ti)+floor(param.seg_len/2))+ (FF.*hann(param.seg_len))';
%% %     ti
%% end
%% % ori_sig = x1(basicTF.fs*time_stamp*(200-1)-100:basicTF.fs*time_stamp*(260-1)+99); ori_sig = ori_sig(201:end-200)';
%% s_factor = param.seg_len/basicTF.hop/2;
%% ma_sig = ma_sig./s_factor; %ma_sig = ma_sig(201:end-200);
%% fe_sig = fe_sig./s_factor; %fe_sig = fe_sig(201:end-200);
%% res = x-ma_sig-fe_sig;
%% % t_index = (1:1000)./100;
