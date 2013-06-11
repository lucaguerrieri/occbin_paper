function chain=markov_luca(P,n,start_pos)

chain = zeros(n,1);
chain(1) = start_pos;

Psum = cumsum(P,2);


for i=2:n
    thisdraw = rand(1);
    chain(i) = min(find(thisdraw<Psum(chain(i-1),:)));    
end

