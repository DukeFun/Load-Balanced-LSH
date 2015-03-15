function T = lshins(T,x,ind)
% T = lshins(T,X)
%
%     insert data (columns of X) into T
% 
% NOTE: LSH will only index the data - not store it! You need to keep
% around the original data, in order to go back from indices to actual
% points, if that's what you want to do.
%
% T = lshins(T,X,IND)
%   instead of assuming that columns of X have indices 1..size(X,2), uses IND
%    
%
% (C) Greg Shakhnarovich, TTI-Chicago (2008)


% fields of T:
% buckets : bukets(j,:) is the hash key of bucket j
% Index : Index{j} contains indices of data in bucket j
% count : count(j) contains the size of bucket j

if (nargin < 3 || isempty(ind))
  ind=1:size(x,2);
end


% insert in each table
for j=1:length(T)   %从1到hash表的个数
  
  % the # of buckets before new data 行
  oldBuckets=size(T(j).buckets,1);
  
  % find, for each data point, the corresp. bucket
  % bucket numbers are represented as arrays of uint8
  % note: 0/1 keys in LSH become 128/129 N*k维hash值矩阵
  buck = findbucket(T(j).type,x,T(j).I);
  
  % now x(:,n) goes to bucket with key uniqBuck(bID(n))
  % 输入数据x(:,n)由uniqBuck(bID(n))表示
  [uniqBuck,ib,bID] = unique(buck,'rows');  % 注意unique函数的返回值
  keys = lshhash(uniqBuck);
  
  if (T(j).verbose > 0)
    fprintf(2,'%d distinct buckets\n',length(ib));
  end
  
  % allocate space for new buckets -- possibly excessive，开辟空间
  %length(ib)代表uniq操作后的数据个数，k代表hash表中的关键字个数
  T(j).buckets=[T(j).buckets; zeros(length(ib),T(j).I.k,'uint8')];
  
  newBuckets=0;
  
  for b=1:length(ib)
    % find which data go to bucket uniqBuck(b)
    thisBucket = find(bID==bID(ib(b)));     %返回与查询值具有相同值的数据的行ID
    
    % find out if this bucket already has anything
    % first, which bucket is it?
    ihash = T(j).bhash{keys(b)}; % possible matching buckets
    if (isempty(ihash)) % nothing matches
      isb = [];
    else % may or may not match
      isb = ihash(find(all(bsxfun(@eq,uniqBuck(b,:),T(j).buckets(ihash,:)),2)));
    end
    
    % note: this search is the most costly operation
    %isb = find(all(bsxfun(@eq,uniqBuck(b,:),T(j).buckets),2));
    
    if (~isempty(isb)) 
      % adding to an existing bucket.
      oldcount=length(T(j).Index{isb}); % # elements in the bucket prior
                                        % to addition
      newIndex = [T(j).Index{isb}  ind(thisBucket)];    %合并 T(j).Index{isb}选中cell中的一个单元
    else
      % creating new bucket
      newBuckets=newBuckets+1;
      oldcount=0;
      isb = oldBuckets+newBuckets;
      T(j).buckets(isb,:)=uniqBuck(b,:);
      T(j).bhash{keys(b)} = [T(j).bhash{keys(b)}; isb];
      newIndex = ind(thisBucket);
    end
    
%     % 当某个bucket超过规定阈值B时的处理过程
%     % if there is a bound on bucket capacity, and the bucket is full,
%     % keep a random subset of B elements (note: we do this rather than
%     % simply skip the new elements since that could introduce bias
%     % towards older elements.)
%     % There is still a bias since older elements have more chances to get
%     % thrown out.
%     if (length(newIndex) > T(j).B)
%       rp=randperm(length(newIndex));
%       newIndex = newIndex(rp(1:T(j).B));
%     end
    % ready to put this into the table
    T(j).Index{isb}= newIndex;
    % update distinct element count
    T(j).count = T(j).count + length(newIndex)-oldcount;
  end
  
   %计算虚拟中心参数
    VC(j).bu = zeros(size(x,1), size(T(j).buckets,1));%虚拟中心
    for c = 1:size(T(j).buckets,1)
		
		VC(j).dist(c) = 0;	%最大距离
		VC(j).newindex(c) = 0;	%
		VC(j).index(c) = 0;	%edge point在数据集中的索引
		
		%if(length(T(j).buckets(c)) >= T(j).B)
		virCenter = zeros(size(x,1), 1);
		inCell = T(j).Index(c);
		idx = inCell{1,1};
		for ci = 1:size(idx,2)
			virCenter = virCenter + x(:, idx(ci));
		end
		virCenter = virCenter / length(idx);
		VC(j).bu(:, c) = virCenter;
        
        if(length(idx) > T(j).B)
            for ci = 1:size(idx,2)
                tempDist = sum( ( x(:, idx(ci) ) - virCenter(:, 1)) .* ( x( :, idx(ci) ) - virCenter(:, 1)), 1);
                if(tempDist >= VC(j).dist(c))
                    VC(j).dist(c) = tempDist;
                    VC(j).newindex(c) = ci;
                    VC(j).index(c) = idx(ci);
                end
            end	
        else
            VC(j).index(c) = 0;
        end
    end
    
    for c = 1:size(T(j).buckets,1)
        
        %新思路，只关心当前bucket
        fprintf(2, 'Processing the %d bucket in the %d Hash Table \n', c, j);
		inCell = T(j).Index(c);
		idx = inCell{1,1};
		while(length(idx) > T(j).B)
            
            VC(j).balance(c) = 0;
            if(VC(j).index(c) == 0)
                for ci = 1:length(idx)
                    tempDist = sum( ( x(:, idx(ci) ) - VC(j).bu(:, c)) .* ( x( :, idx(ci) ) - VC(j).bu(:, c)), 1);
                    if(tempDist >= VC(j).balance(c))
                        VC(j).balance(c) = tempDist;
                        VC(j).newindex(c) = ci;
                        VC(j).index(c) = idx(ci);
                    end
                end
            end
            
            if( c == size(T(j).buckets, 1))
                inCellPlus = T(j).Index(1);
                idxPlus = inCellPlus{1,1};
                inCellPlus{1,1} = [idxPlus VC(j).index(c)];
                T(j).Index(1) = inCellPlus;
            else
                inCellPlus = T(j).Index(c + 1);
                idxPlus = inCellPlus{1,1};
                inCellPlus{1,1} = [idxPlus VC(j).index(c)];
                T(j).Index(c + 1) = inCellPlus;
            end
            
            idx( VC(j).newindex(c) ) = [];
            inCell{1,1} = idx;
			%inCell{1,1} = [ idx(1 : VC(j).newindex - 1) idx(VC(j).newindex + 1 : length(idx)) ];
			T(j).Index(c) = inCell;
			
			inCell = T(j).Index(c);
			idx = inCell{1,1};
            VC(j).balance(c) = 0;
			for ci = 1:length(idx)
				tempDist = sum( ( x(:, idx(ci) ) - VC(j).bu(:, c)) .* ( x( :, idx(ci) ) - VC(j).bu(:, c)), 1);
				if(tempDist >= VC(j).balance(c))
					VC(j).balance(c) = tempDist;
					VC(j).newindex(c) = ci;
					VC(j).index(c) = idx(ci);
				end
			end
			
			if( c == size(T(j).buckets, 1))
				c = 0;
            end
			
			if( c == 0)
				c = c + 1;
			end
        end
    end
    
  % we may not have used all of the allocated bucket space
  T(j).buckets=T(j).buckets(1:(oldBuckets+newBuckets),:);
  if (T(j).verbose > 0)
    fprintf(2,'Table %d adding %d buckets (now %d)\n',j,newBuckets,size(T(j).buckets,1));
    fprintf(2,'Table %d: %d elements\n',j,T(j).count);
  end
end



