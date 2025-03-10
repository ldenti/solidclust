CPU_CC        := $(CC)
CPU_OPTS      := --std=c89 -Wall -Werror -g #-O3
CPU_MODULES   := exe lib bundled/klib bundled/xxhash bundled/xoshiro
CPU_LIBS	  := -lm -lz # -lpthread
CPU_BUILD_DIR := build
CPU_OBJ_DIR   := $(CPU_BUILD_DIR)/objs
CPU_SRC_DIR   := $(foreach item,$(CPU_MODULES),$(item)/src)
CPU_INC_DIR   := $(foreach item,$(CPU_MODULES),$(item)/include)
CPU_SRC       := $(foreach sdir,$(CPU_SRC_DIR),$(wildcard $(sdir)/*.c))
CPU_INC       := $(addprefix -I,$(CPU_INC_DIR))
CPU_BASE      := $(foreach sdir,$(CPU_SRC_DIR),$(notdir $(wildcard $(sdir)/*.c)))
CPU_OBJ       := $(patsubst %.c,$(CPU_OBJ_DIR)/%.o,$(CPU_BASE))

.PHONY: all checkdirs clean

all: checkdirs build/isonclust

build/isonclust: $(CPU_OBJ)
	$(CPU_CC) $^ -o $@ $(CPU_LIBS)

checkdirs: $(CPU_OBJ_DIR)

$(CPU_OBJ_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(CPU_BUILD_DIR) 

.SECONDEXPANSION:
PERCENT = %

$(CPU_OBJ): %.o : $$(filter $$(PERCENT)/$$(notdir %).c, $(CPU_SRC))
	$(CPU_CC) $(CPU_OPTS) $(CPU_INC) -c $< -o $@

$(DPU_OBJ): %.o : $$(filter $$(PERCENT)/$$(notdir %).c, $(DPU_SRC))
	$(DPU_CC) -DNR_TASKLETS=$(NR_TASKLETS) -DSTACK_SIZE_DEFAULT=$(STACK_SIZE_DEFAULT) $(DPU_OPTS) $(DPU_INC) -c $< -o $@
