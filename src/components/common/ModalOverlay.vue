<!-- 通用全屏弹层基座：Teleport + 滚动锁定 + Esc 全局关闭 + backdrop 点击关闭 + 焦点还原 -->
<script setup lang="ts">
import { onBeforeUnmount, onMounted } from 'vue'

import { lockScrollOverflow, unlockScrollOverflow } from '@/utils/scroll'

const emit = defineEmits<{ close: [] }>()

let previousFocus: HTMLElement | null = null

const onKeydown = (e: KeyboardEvent) => {
  if (e.key === 'Escape') emit('close')
}

onMounted(() => {
  previousFocus = document.activeElement as HTMLElement | null
  lockScrollOverflow()
  window.addEventListener('keydown', onKeydown)
})

onBeforeUnmount(() => {
  unlockScrollOverflow()
  window.removeEventListener('keydown', onKeydown)
  previousFocus?.focus?.()
})
</script>

<template>
  <Teleport to="body">
    <div class="modal-overlay" @click.self="emit('close')">
      <slot />
    </div>
  </Teleport>
</template>

<style scoped>
.modal-overlay {
  position: fixed;
  inset: 0;
  z-index: 1080;
  background: rgba(0, 0, 0, 0.42);
  display: flex;
  align-items: flex-start;
  justify-content: center;
  padding-top: 14vh;
}
</style>
