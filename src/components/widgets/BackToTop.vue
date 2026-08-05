<!-- BackToTop.vue -->
<template>
  <button
    v-show="showBackToTop"
    class="back-to-top d-flex align-items-center justify-content-center"
    @click="handleClick"
    :aria-label="t('backToTop')"
    @touchstart.prevent.stop="handleTouchStart"
    @touchmove.prevent.stop="handleTouchMove"
    @touchend.prevent.stop="handleTouchEnd"
    :style="{ top: buttonTop + 'px' }"
  >
    <i class="fas fa-arrow-up"></i>
  </button>
</template>

<script setup lang="ts">
import { onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'
import { useFloatingButton } from '@/composables/useFloatingButton'

const { t } = useI18n()

const showBackToTop = ref(false)
const scrollObserver = ref<IntersectionObserver | null>(null)
const scrollSentinel = ref<HTMLDivElement | null>(null)

function backToTop() {
  const reduceMotion =
    typeof window.matchMedia === 'function' &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches
  window.scrollTo({ top: 0, behavior: reduceMotion ? 'auto' : 'smooth' })
}

const {
  isDragging,
  buttonTop,
  dispatchBaseTop,
  onTouchStart,
  onTouchMove,
  onTouchEnd,
  subscribe,
  unsubscribe,
} = useFloatingButton({
  sourceId: 'btt',
  defaultTop: (typeof window !== 'undefined' ? window.innerHeight : 1024) - 100,
  mode: 'match',
  onRelease: backToTop,
})

const handleClick = () => {
  if (!isDragging.value) backToTop()
}
const handleTouchStart = onTouchStart
const handleTouchMove = onTouchMove
const handleTouchEnd = onTouchEnd

function setupScrollObserver() {
  if (scrollObserver.value) return
  const sentinel = document.createElement('div')
  sentinel.style.cssText =
    'position:absolute;left:0;top:0;width:1px;height:181px;pointer-events:none;visibility:hidden'
  document.body.appendChild(sentinel)
  scrollSentinel.value = sentinel

  scrollObserver.value = new IntersectionObserver(
    ([entry]) => {
      showBackToTop.value = !entry.isIntersecting
      if (showBackToTop.value) dispatchBaseTop()
    },
    { threshold: 0 },
  )
  scrollObserver.value.observe(sentinel)
}

onMounted(() => {
  setupScrollObserver()
  subscribe()
  buttonTop.value = window.innerHeight - 100
  dispatchBaseTop()
})

onBeforeUnmount(() => {
  if (scrollObserver.value) {
    scrollObserver.value.disconnect()
    scrollObserver.value = null
  }
  if (scrollSentinel.value) {
    scrollSentinel.value.remove()
    scrollSentinel.value = null
  }
  unsubscribe()
})
</script>

<style scoped>
.back-to-top {
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

.back-to-top:hover {
  background-color: var(--primary);
  color: var(--on-primary);
  border-color: transparent;
  box-shadow: var(--shadow-lift);
}

.back-to-top:active {
  transform: scale(0.93);
}
</style>
