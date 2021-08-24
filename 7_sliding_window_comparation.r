
seq_string <- DNAStringSet(my_fasta_string)


consensusMatrix <- consensusMatrix(seq_string)
pwm_0 <- consensusMatrix[DNA_BASES, ]
pwm_1 <- prop.table(pwm_0, 2)


pwm <- makePWM(pwm_1)

slotNames(pwm)
pwm(pwm)
ic(pwm)
consensus(pwm)

bit_score <- ic(pwm)



begin_site<- 1:65
sum_min=10000000
bgein_min=0
for (i in begin_site) {
  begin=i
  begin_2=begin+10
  sum_i = sum(bit_score[begin:begin_2])
  
  if(sum_i < 18){
    sum_min = sum_min
    bgein_min=bgein_min
  } else {
    print(begin)
    begin_base_site = begin
    break
  }
}


finish_site <- length(bit_score):(length(bit_score)-30)
sum_min=10000000
finish_min=0

for (i in finish_site) {
  finish=i
  finish_2=finish-10
  sum_i = sum(bit_score[finish:finish_2])
  
  if(sum_i < 18){
    sum_min = sum_min
    bgein_min=bgein_min
  } else {
    print(finish)
    finish_base_site = finish
    break
  }
}






short_base_length=230
seq_length=finish_base_site-begin_base_site

base_site <- begin_base_site:(seq_length-short_base_length)
sum_min=10000000
start_min=0
end_min=0
for (i in base_site) {
  start=i
  end=start+short_base_length
  sum_i = sum(bit_score[start:end])
  
  if(sum_i < sum_min){
    sum_min = sum_i
    start_min=start
    end_min=end
  } else {
    sum_min = sum_min
    start_min=start_min
    end_min=end_min
    }
}


print(sum_min)
print(start_min)
print(end_min)

cut_start_1=start_min
cut_end_1=start_min+((end_min-start_min)/2)
cut_start_2=start_min+((end_min-start_min)/2)
cut_end_2=end_min
