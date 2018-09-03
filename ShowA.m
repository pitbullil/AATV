function ShowA( A , i )
if nargin==1
    i=121;
end
figure(i);
subplot 221; imagesc(A{1,1}); title('A{1,1}');
subplot 222; imagesc(A{2,1}); title('A{2,1}');
subplot 223; imagesc(A{1,2}); title('A{1,2}');
subplot 224; imagesc(A{2,2}); title('A{2,2}');
end

