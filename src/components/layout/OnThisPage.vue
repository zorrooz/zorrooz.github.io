<template>
  <nav class="on-this-page">
    <div class="otp-header">
      <span class="otp-title">{{ t('tableOfContents') }}</span>
    </div>

    <ul class="otp-list">
      <li
        v-for="(item, idx) in toc"
        :key="idx"
        :class="['otp-item', { active: activeId === item.id }]"
      >
        <a
          class="otp-link"
          role="button"
          tabindex="0"
          @click.prevent="scrollToId(item.id)"
          @keydown.enter.prevent="scrollToId(item.id)"
        >
          <span class="otp-text">{{ item.text }}</span>
        </a>
        <ul v-if="item.children && item.children.length" class="otp-sublist">
          <li
            v-for="(sub, sIdx) in item.children"
            :key="sIdx"
            :class="['otp-subitem', { active: activeId === sub.id }]"
          >
            <a
              class="otp-sublink"
              role="button"
              tabindex="0"
              @click.prevent="scrollToId(sub.id)"
              @keydown.enter.prevent="scrollToId(sub.id)"
            >
              <span class="otp-subtext">{{ sub.text }}</span>
            </a>
          </li>
        </ul>
      </li>
    </ul>
  </nav>
</template>

<script setup lang="ts">
import { nextTick, onBeforeUnmount, onMounted, ref } from 'vue'
import { useI18n } from 'vue-i18n'

interface TocNode {
  id: string
  text: string
  level: number
  children: TocNode[]
}

const props = withDefaults(
  defineProps<{
    containerSelector?: string
    levels?: number[]
    offset?: number
  }>(),
  {
    containerSelector: '.markdown-body',
    levels: () => [2, 3, 4, 5, 6],
    offset: 8,
  },
)

const emit = defineEmits(['navigate'])

const { t } = useI18n()

const toc = ref<TocNode[]>([])
const activeId = ref('')
const otpObserver = ref<MutationObserver | null>(null)
const otpObserverTimer = ref<number | null>(null)
const otpPoller = ref<number | null>(null)
const spyObserver = ref<IntersectionObserver | null>(null)

function cleanupObservers() {
  if (otpObserver.value) {
    otpObserver.value.disconnect()
    otpObserver.value = null
  }
  if (otpObserverTimer.value) {
    clearTimeout(otpObserverTimer.value)
    otpObserverTimer.value = null
  }
  if (otpPoller.value) {
    clearInterval(otpPoller.value)
    otpPoller.value = null
  }
  if (spyObserver.value) {
    spyObserver.value.disconnect()
    spyObserver.value = null
  }
}

function resetToc() {
  toc.value = []
  activeId.value = ''
  cleanupObservers()
  nextTick(() => setupContainerObserver())
}

function refreshToc() {
  buildToc()
  setupScrollSpy()
}

function setupContainerObserver() {
  setTimeout(() => refreshToc(), 100)
  const checkContainer = () => {
    const root = document.querySelector(props.containerSelector)
    if (!root) return
    if (otpPoller.value) {
      clearInterval(otpPoller.value)
      otpPoller.value = null
    }
    if (!otpObserver.value) {
      otpObserver.value = new MutationObserver(() => {
        clearTimeout(otpObserverTimer.value ?? undefined)
        otpObserverTimer.value = window.setTimeout(() => refreshToc(), 100)
      })
      otpObserver.value.observe(root, {
        childList: true,
        subtree: true,
        attributes: true,
        attributeFilter: ['id'],
      })
    }
    refreshToc()
  }
  checkContainer()
  if (!otpPoller.value && !document.querySelector(props.containerSelector))
    otpPoller.value = window.setInterval(checkContainer, 200)
}

function getHeadingText(h: Element) {
  try {
    const clone = h.cloneNode(true) as Element
    clone.querySelectorAll('.heading-anchor')?.forEach((a) => a.remove())
    return (clone.textContent || '').replace(/\s*#\s*$/, '').trim()
  } catch {
    return (h.textContent || '').replace(/\s*#\s*$/, '').trim()
  }
}

function buildToc() {
  const root = document.querySelector(props.containerSelector)
  if (!root) {
    toc.value = []
    return
  }

  const selector = props.levels.map((l) => `h${l}`).join(',')
  const headings = Array.from(root.querySelectorAll(selector))

  if (headings.length === 0) {
    toc.value = []
    return
  }

  headings.forEach((h) => {
    if (!h.id) {
      let safeId = h.textContent
        .trim()
        .toLowerCase()
        .replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g, '')
        .replace(/\s+/g, '-')

      if (!safeId) {
        safeId = `section-${Math.random().toString(36).substring(2, 9)}`
      }

      let finalId = safeId
      let count = 1
      while (document.getElementById(finalId)) {
        finalId = `${safeId}-${count++}`
      }
      h.id = finalId
    }
  })

  const levelSet = new Set(props.levels)
  const topLevel = Math.min(...props.levels)
  const secondLevel = topLevel + 1
  const tocList: TocNode[] = []
  let currentTop: TocNode | null = null
  for (const h of headings) {
    const level = parseInt(h.tagName.substring(1), 10)
    if (!levelSet.has(level)) continue
    const node: TocNode = { id: h.id, text: getHeadingText(h), level, children: [] }
    if (level === topLevel) {
      tocList.push(node)
      currentTop = node
    } else if (currentTop && level >= secondLevel) currentTop.children.push(node)
    else tocList.push(node)
  }
  toc.value = tocList
}

function setupScrollSpy() {
  if (spyObserver.value) {
    spyObserver.value.disconnect()
    spyObserver.value = null
  }
  const root = document.querySelector(props.containerSelector)
  if (!root) return

  const selector = props.levels.map((l) => `h${l}`).join(',')
  const headings = Array.from(root.querySelectorAll(selector))
  if (headings.length === 0) return

  spyObserver.value = new IntersectionObserver(
    (entries) => {
      for (const entry of entries) {
        if (entry.isIntersecting) activeId.value = (entry.target as HTMLElement).id
      }
    },
    { rootMargin: `-${props.offset}px 0px -60% 0px`, threshold: 0 },
  )
  headings.forEach((h) => spyObserver.value?.observe(h))
}

function scrollToId(id: string) {
  emit('navigate', id)

  const el = document.getElementById(id)
  if (!el) {
    return
  }
  const top = el.getBoundingClientRect().top + window.scrollY - props.offset
  const reduceMotion =
    typeof window.matchMedia === 'function' &&
    window.matchMedia('(prefers-reduced-motion: reduce)').matches
  const doScroll = () => {
    window.scrollTo({
      top,
      behavior: reduceMotion ? 'auto' : 'smooth',
    })
  }

  try {
    const bodyOverflow = document.body && document.body.style && document.body.style.overflow
    if (bodyOverflow === 'hidden') {
      setTimeout(doScroll, 80)
    } else {
      doScroll()
    }
  } catch {
    doScroll()
  }
}

onMounted(() => {
  buildToc()
  setupScrollSpy()

  nextTick(() => {
    setupContainerObserver()
  })
})

onBeforeUnmount(() => {
  cleanupObservers()
})

defineExpose({ refreshToc, resetToc })
</script>

<style scoped>
.on-this-page {
  padding: 0.25rem 0;
  font-size: var(--text-base);
}

.otp-header {
  display: flex;
  align-items: center;
  gap: var(--sp-2);
  margin-bottom: 0.75rem;
}

.otp-title {
  font-weight: 600;
  font-size: 12px;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: var(--fg-3);
}

.otp-list {
  list-style: none;
  padding: 0;
  margin: 0;
}

.otp-item {
  margin: 1px 0;
  padding-left: 0;
}

.otp-link {
  display: flex;
  align-items: center;
  gap: var(--sp-2);
  padding: 6px 10px;
  color: var(--fg-2);
  text-decoration: none;
  font-size: var(--text-sm);
  border-radius: var(--radius-sm);
  transition:
    color 0.14s ease,
    background-color 0.14s ease;
  position: relative;
  line-height: 1.5;
  cursor: pointer;
}

.otp-link:hover {
  color: var(--primary);
  background-color: var(--surface-2);
}

.otp-item.active > .otp-link {
  color: var(--primary);
  font-weight: 600;
}

.otp-item.active > .otp-link::before {
  content: '';
  position: absolute;
  left: 0;
  top: 50%;
  transform: translateY(-50%);
  width: 2px;
  height: 12px;
  border-radius: 99px;
  background: var(--primary);
}

.otp-sublist {
  list-style: none;
  padding-left: 12px;
  margin: 2px 0;
  border-left: 1px solid var(--line);
}

.otp-subitem .otp-sublink {
  display: block;
  padding: 4px 10px;
  color: var(--fg-3);
  text-decoration: none;
  font-size: var(--text-sm);
  border-radius: var(--radius-sm);
  transition:
    color 0.14s ease,
    background-color 0.14s ease;
  cursor: pointer;
  line-height: 1.5;
}

.otp-sublink:hover {
  color: var(--primary);
  background-color: var(--surface-2);
}

.otp-subitem.active > .otp-sublink {
  color: var(--primary);
  font-weight: 600;
}
</style>
