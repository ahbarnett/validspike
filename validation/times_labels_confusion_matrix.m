function [confusion_matrix T1 T2 T2w permL2]=times_labels_confusion_matrix(T1,L1,T2,L2,opts)
% TIMES_LABELS_CONFUSION_MATRIX - confusion matrix between two {t,l} event lists
%
% [P t1 t2 t2w permL2] = times_labels_confusion_matrix(T1,L1,T2,L2,opts)
%  returns confusion matrix P of size (K1+1)-by-(K2+1) where Ki is the number
%  of types in the label list Li, i=1,2. The last row and column of P count
%  the unclassified or missed events in one of the event lists. The list L2 is
%  permuted in the best fashion to match.
%
% Inputs:
%  T1,L1 - list of times (real) and labels (in 1..K1) for firing event list 1
%  T2,L2 - list of times (real) and labels (in 1..K2) for firing event list 2
%  opts  - (optional) controls various parameters such as:
%          opts.max_matching_offset : max time difference to be considered a
%                                     match (positive real, in whatever units
%                                     T1,T2 are in. Default 3.0)
%          opts.internal_run_code (default='main') -- an internal parameter 
% Outputs:
%  P - confusion matrix (extended, ie size (K1+1)-by-(K2+1)), best permuted.
%  t1,t2,t2w - the lists of unmatched times in T1, unmatched times in T2, and
%            wrongly matched times in T2.
%  permL2 - best permutation of L2, ie the labels permL2(L2) are best matching
%           against L1, given the times.
%
% Notes: has 3 passes, in order to find best permutation of labels.
%  Time offsets are only compared to the nearest integer in whatever units T1
%  and T2 are in, due to the integer stepping in offset.
% 
% Run without arguments gives a self-test

% JFM. doc & extra outputs & self-test by AHB 5/15/15
% ahb: if the matching is incorrect, but a correct matching exists within
%  the allowed max_matching_offset, it will not be found. It is "greedy" wrt
%  offset times. -- jfm put in an attempted fix on 5/21/15.
% bug fixes by jfm 5/21/2015
% 5/22/2015 JFM added a crazy third pass because we first need to permute
% the labels for a best match, then run it with two passes as described
% below, and then finally undo the permutation.
% 6/12/15: ahb switched to perm of L2 not L1, and permuted C-matrix output.
% 7/20/15: ahb sped up the offending O(N1.N2.max_matching_offset) lines,
%  pinned down issue of mismatches due to greedy-in-time-offset.
% 8/14/15: ahb switched to use bestcolpermconfmat

% todo: check "unexpected problem" when there's no labels of a certain type
%       (see testcase below).

% run the test if no arguments
if nargin==0, test_times_labels_confusion_matrix; return; end

%default parameters
if nargin<5, opts = []; end
if (~isfield(opts,'max_matching_offset')) opts.max_matching_offset=3; end;
if (~isfield(opts,'internal_run_code')) opts.internal_run_code='main'; end;

% handle some trivial cases to prevent later crashes...
if isempty(T1), confusion_matrix = histc(L2,1:max(L2));
  confusion_matrix = [confusion_matrix(:)', 0];        % row. T2 all unmatched
  T2w = []; permL2 = 1:max(L2); return
end
if isempty(T2), confusion_matrix = histc(L1,1:max(L1));
  confusion_matrix = [confusion_matrix(:); 0];        % col
  T2w = []; permL2 = 1:max(L2); return
end
 

[T1,i] = sort(T1); L1 = L1(i); [T2,i] = sort(T2); L2 = L2(i); % time sort both lists

if (strcmp(opts.internal_run_code,'main'))
	%first we need to find the best permutation as described a bit above
	opts2=opts; opts2.internal_run_code='skip_first_pass';
	CCC=times_labels_confusion_matrix(T1,L1,T2,L2,opts2);
	CCC=CCC(1:end-1,1:end-1); %exclude unclassified row/column
        [permL2] = bestcolpermconfmat(CCC);
	[K1,K2]=size(CCC);
	opts2=opts; opts2.internal_run_code='normal'; % normal means don't find the permutation again
        ii = (L2>=1 & L2<=K2);                % classified indices in L2
        L2(ii) = permL2(L2(ii));              % reshuffle L2's classified labels
	[DDD,T1,T2,T2w]=times_labels_confusion_matrix(T1,L1,T2,L2,opts2);
	if ((size(DDD,1)~=K1+1)||(size(DDD,2)~=K2+1))
		disp('sizes of two confusion matrices:');
                disp (size(DDD)); disp ([K1+1,K2+1]);
		error('Unexpected problem: confusion matrix size problem - maybe not all labels are represented? Can fix');
	end;
	confusion_matrix=DDD; % ahb removed permuting conf-mat since it's already permuted! :
	return;
end;

skip_first_pass=strcmp(opts.internal_run_code,'skip_first_pass');

% the number of labels is the maximum of L1 and L2, respectively
num_labels1=max(L1);
num_labels2=max(L2);
T2w = [];   % AHB keep list of wrongly-matched

%initialize the output confusion matrix
confusion_matrix=zeros(num_labels1+1,num_labels2+1); %the last row and column are for unclassified events
	
%first we handle the pairs that are 0 close, then 1 close, then 2 close,
%all the way up to max_offset close.
%5/21/15: Except we want to first give the benefit of the doubt and match all events
%that agree (in the first pass). This is to handle the situation pointed out by ahb where
%events are close enough to match, but they happen to be offset a bit so
%they get paired up with the wrong events (jfm - 5/12/2015)
for pass=1:2 %first pass handles case where matching events should be given priority
	if ((pass~=1)||(~skip_first_pass))            % three passes in total
		for offset=0:opts.max_matching_offset
		  inds1_to_remove=zeros(1,length(T1)); %once paired up we remove those events
		  inds2_to_remove=zeros(1,length(T2)); %once paired up we remove those events
		  if (length(T2)>0) %jfm 5/27/15
                          % precompute cond_for_pass_1 for L1=1...K1 here?
                          ptr2 = 1;  % moving index to sorted T2 list
			  for it=1:length(T1) %this loop makes things slow (should we fix?)
				%find the indices of T2/L2 which should be paired with the event at T1(it)
                                [ii,ptr2] = indexlist(T2,T1(it),offset,ptr2);
                                %ii = 1:numel(T2); % uncomment for original slow method
				if (pass==1)
                                  condition_for_pass_1=(L2(ii)==L1(it)); %first pass handles case where matching events should be given priority
                                                 % could be sped up by precomputing for L1 = 1..K1
                                else, condition_for_pass_1=ones(size(ii)); %fixed by jfm on 5/27/15, moved ahb 7/20/15
				end
				inds=ii(find((abs(T2(ii)-T1(it))<=offset).*(inds2_to_remove(ii)==0).*condition_for_pass_1)); % ahb sped up
				if (length(inds)>0)
				  ind0=inds(1);  % Handles case of a double match
				  confusion_matrix(L1(it),L2(ind0))=confusion_matrix(L1(it),L2(ind0))+1; %increment the entry in the confusion matrix
				  inds1_to_remove(it)=1; %we've handled the event, so let's remove it!
				  inds2_to_remove(ind0)=1; %we've handled the event, so let's remove it!
				  if L2(ind0)~=L1(it), T2w = [T2w T1(it)]; end % AHB - keep track of the wrongly classified events
				end;
			  end;
			  %Now let's remove the events that were marked above
			  T1=T1(inds1_to_remove==0); 
			  L1=L1(inds1_to_remove==0);
			  T2=T2(inds2_to_remove==0);
			  L2=L2(inds2_to_remove==0);
		  end;
		end;
	end;
end;
%The remaining are unclassified
for it=1:length(T1)
  confusion_matrix(L1(it),num_labels2+1)=confusion_matrix(L1(it),num_labels2+1)+1;
end
for it=1:length(T2)
  confusion_matrix(num_labels1+1,L2(it))=confusion_matrix(num_labels1+1,L2(it))+1;
end

%%%%%%%
function [ii,ptr2] = indexlist(t2,t1,off,ptr2) % find indices of t2 for which
% abs(t2-t1)<off, where t1 is a scalar. Update ptr2 which gives the rough
% center value of this index. Barnett 7/20/15
N = numel(t2);
i = ptr2; while(t2(i)>t1-off & i>1), i=i-1; end  % go down until t2's too early
j = ptr2; while(t2(j)<t1+off & j<N), j=j+1; end  % go up until t2's too late
ii = i:j;
ptr2 = round((i+j)/2);


%%%%%%%%%%%%%%%
function test_times_labels_confusion_matrix
% Barnett 5/15/15. Better tests 7/20/15

if 1
  times_labels_confusion_matrix([10 20],[1 2], [10 20 30 40], [1 2 3 4])
  
  T=[10 20 30]; L=[1 2 3]; S=[10 20]; M=[1 2];  % debug of 5/27/15
  [P t1 t2 t2w] = times_labels_confusion_matrix(T,L,S,M)
  [P t1 t2 t2w] = times_labels_confusion_matrix(S,M,T,L)
  [P t1 t2 t2w] = times_labels_confusion_matrix(T,L,[10 21],M)
end

if 1
T = [10 20 30]; L = [1 2 2];  % 3 spikes, times far from each other

'should show t2w only:'
L2 = [1 2 3]; i = randperm(3); T2 = T(i); L2 = L2(i);
[P t1 t2 t2w] = times_labels_confusion_matrix(T,L,T2,L2)

'should show t1 only:'
T2 = T(1:2); L2 = [1 2];
[P t1 t2 t2w] = times_labels_confusion_matrix(T,L,T2,L2)

'should show t1 and t2, but no t2w:'
T2 = T; T(3) = 40; L2 = L;
[P t1 t2 t2w] = times_labels_confusion_matrix(T,L,T2,L2)

'extra:'
[P t1 t2 t2w]=times_labels_confusion_matrix([10 20],[1 2],[10 20 30],[1 2 3])
'missing:'
[P t1 t2 t2w]=times_labels_confusion_matrix([10 20 30],[1 2 3],[10 20],[1 2])
  
% two nearby times, consistent, but it matches wrong and says they are
% wrong: (should be fixed -- jfm)
[P t1 t2 t2w] = times_labels_confusion_matrix([10 12],[1 2],[11 10],[1 2]);
end

if 0  % jfm's debugging
[P t1 t2 t2w] = times_labels_confusion_matrix([10 12 15 20],[1 1 2 2],[11 14 16 21],[2 2 1 1],struct('internal_run_code','main'));
disp(P);
T1=[10 1 2 12 15 18 80 80 80 21 52 54]; L1=[1 1 1 1 2 2 2 2 2 3 3 3];
T2=[11 1 2 14 16 21 80 80 80 24 52 54]; L2=[2 2 2 2 3 3 3 3 3 1 1 1];
[P t1 t2 t2w] = times_labels_confusion_matrix(T1,L1,T2,L2);
disp(P);
[P t1 t2 t2w] = times_labels_confusion_matrix(T1,L1,T2,L2,struct('internal_run_code','normal')); %skip the permutation step and it should give wrong answer
disp(P);
end

if 1 % these caused it to crash (should be fixed -- jfm)
[P t1 t2 t2w] = times_labels_confusion_matrix([10 20],[1 1],[12],[1])
[P t1 t2 t2w] = times_labels_confusion_matrix([10 20 30],[1 1 2],[12 20],[1 2])
end

if 1 % "unexpected problem" testcase
T = [2905.3       2908.5       2909.6]; L=[2 1 2];   %L = [7     6     7];
T2 = [2907.4       2908.1       2911.5]; L2=[1 2 2]; %L2=[6     2     2];
T = [T 1e4+(10:10:100)]; T2 = [T2 1e4+(10:10:100)]; L=[L ones(1,10)]; L2=[L2 ones(1,10)];
o.max_matching_offset = 10;
[Q tmiss tfals twrng Pfound] = times_labels_confusion_matrix(T,L,T2,L2,o)
end

if 1 % Alex's test on random long list of times & labels...
N=1e4; K=7; T = 10*N*rand(1,N); L = randi(K,1,N); % note firing rate matters
P = randperm(K), L2 = P(L); T2 = T + 3.0*(2*rand(1,N)-1);  % permute & jitter
o.max_matching_offset = 5;
tic; [Q tmiss tfals twrng Pfound] = times_labels_confusion_matrix(T,L,T2,L2,o);
toc
Q
fprintf('confusion matrix error frac = %.3g\n',(N-sum(diag(Q)))/N)
% note that greedy matching as increase t offset can lead to non-exhaustive
% searching, hence errors here, even though there exists a complete match!
fprintf('best permutation error frac = %.3g\n',sum(Pfound(P)~=1:K)/K)
end
%keyboard
