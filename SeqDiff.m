function [features] = SeqDiff(sequence,idealSeq,partition)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%sequence='CGATCGATCGAATACCGAATGACAGT';
%partition=2;
%idealSeq='GGGG';

w=length(sequence);
distances=zeros(w,1);
idealSeq = convertStringsToChars(idealSeq(1));
idealSeqLength=length(idealSeq);
%%LBP
%scale = 2.^[7 6 5; 0 -inf 4; 1 2 3];
%iterate whole sequence
for i=1:w
    decNum=0;
    t=0;
    %match full template
    for j=1:idealSeqLength
        %calculate similarity for right side sequences
        if (i+idealSeqLength-j)<=w
            idealSeqBase=idealSeq(idealSeqLength+1-j);
            seqBase=sequence(i+idealSeqLength-j);
            if(idealSeqBase=='A')
                %AA=11, AT=00, AC=01, AG=10
                if(seqBase=='A')
                    decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                elseif(seqBase=='C')
                    decNum=decNum+1*2^t;
                    %decNum=decNum+1*2^(t+1);
                elseif(seqBase=='G')
                    %decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                end
            elseif(idealSeqBase=='T')
                %TA=00, TT=11, TC=01, TG=10
                if(seqBase=='T')
                    decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                elseif(seqBase=='C')
                    decNum=decNum+1*2^t;
                    %decNum=decNum+1*2^(t+1);
                elseif(seqBase=='G')
                    %decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                end
            elseif(idealSeqBase=='C')
                %CA=00, CT=01, CC=11, CG=10
                if(seqBase=='T')
                    decNum=decNum+1*2^t;
                    %decNum=decNum+1*2^(t+1);
                elseif(seqBase=='C')
                    decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                elseif(seqBase=='G')
                    %decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                end
            elseif(idealSeqBase=='G')
                %GA=00, GT=01, GC=10, GG=11
                if(seqBase=='T')
                    decNum=decNum+1*2^t;
                    %decNum=decNum+1*2^(t+1);
                elseif(seqBase=='C')
                    %decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                elseif(seqBase=='G')
                    decNum=decNum+1*2^t;
                    decNum=decNum+1*2^(t+1);
                end
            end
            t=t+2;
        end
    end
    distances(i)=decNum;
end

distances2D=0;

step=floor(w/partition);
binSize=power(2,idealSeqLength*2);
binc = 0:binSize-1;
for i=1:partition
    startPos=(i-1)*step+1;
    endPos=i*step;
    partitionHist=hist(distances(startPos:endPos),binc);
    distances2D(i,1:binSize)=partitionHist;
end
features = distances2D(:);
end



