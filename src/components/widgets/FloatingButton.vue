<template>
  <button
    v-show="show"
    class="floating-btn d-flex align-items-center justify-content-center"
    @click="handleClick"
    :aria-label="ariaLabel"
    @touchstart.prevent.stop="onTouchStart"
    @touchmove.prevent.stop="onTouchMove"
    @touchend.prevent.stop="onTouchEnd"
    :style="{ top: buttonTop + 'px' }"
  >
    <slot></slot>
  </button>
</template>

<script setup lang="ts">
import { onBeforeUnmount, onMounted, watch } from 'vue'
import { useFloatingButton } from '@/composables/useFloatingButton'

const props = withDefaults(
  defineProps<{
    sourceId: string
    defaultTop?: number
    mode?: 'match' | 'stack'
    onRelease?: () => void
    ariaLabel: string
    show?: boolean
  }>(),
  {
    mode: 'match',
    show: true,
  },
)

const release = () => {
  props.onRelease?.()
}

const {
  isDragging,
  buttonTop,
  clampTop,
  dispatchBaseTop,
  onTouchStart,
  onTouchMove,
  onTouchEnd,
  subscribe,
  unsubscribe,
} = useFloatingButton({
  sourceId: props.sourceId,
  defaultTop:
    props.defaultTop ?? (typeof window !== 'undefined' ? window.innerHeight : 1024) - 100,
  mode: props.mode,
  onRelease: release,
})

const handleClick = () => {
  if (!isDragging.value) release()
}

function onResize() {
  buttonTop.value = clampTop(buttonTop.value)
}

watch(
  () => props.show,
  (visible) => {
    if (visible) dispatchBaseTop()
  },
)

onMounted(() => {
  window.addEventListener('resize', onResize, { passive: true })
  subscribe()
  dispatchBaseTop()
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', onResize)
  unsubscribe()
})
</script>

<style scoped>
.floating-btn {
  position: fixed;
  right: 28px;
  width: 40px;
  height: 40px;
  background-color: var(--surface);
  color: var(--fg-2);
  border: 1px solid var(--line);
  border-radius: 50%;
  font-size: 13px;
  cursor: pointer;
  box-shadow: var(--shadow-soft);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
  transition:
    background-color 0.14s ease,
    box-shadow 0.14s ease,
    transform 0.14s ease,
    color 0.14s ease,
    border-color 0.14s ease;
}

.floating-btn:hover {
  background-color: var(--primary);
  color: var(--on-primary);
  border-color: transparent;
  box-shadow: var(--shadow-lift);
}

.floating-btn:active {
  transform: scale(0.93);
}
</style>
