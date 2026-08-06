<!-- AppFooter.vue -->
<script setup lang="ts">
import { useI18n } from 'vue-i18n'
import { SITE } from '@/config'

const { t } = useI18n()

const currentYear = new Date().getFullYear()
const copyrightYears =
  SITE.startYear && SITE.startYear < currentYear
    ? `${SITE.startYear} - ${currentYear}`
    : `${currentYear}`
</script>
<template>
  <footer class="site-footer">
    <div class="container">
      <div class="site-footer__inner">
        <p class="footer-copy">&copy; {{ copyrightYears }} {{ SITE.author }}</p>
        <div class="footer-designed">
          <p class="footer-designed__text">
            {{ t('designedByPrefix')
            }}<strong class="footer-designed__name">{{ SITE.author }}</strong
            >{{ t('designedBySuffix') }}
          </p>
          <div class="footer-designed__icons">
            <a
              :href="SITE.github"
              target="_blank"
              rel="noopener noreferrer"
              class="footer-link"
              :aria-label="'GitHub'"
            >
              <svg
                class="footer-icon"
                width="16"
                height="16"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                stroke-width="1.75"
                stroke-linecap="round"
                stroke-linejoin="round"
                aria-hidden="true"
              >
                <path
                  d="M9 19c-5 1.5-5-2.5-7-3m14 6v-3.87a3.37 3.37 0 0 0-.94-2.61c3.14-.35 6.44-1.54 6.44-7A5.44 5.44 0 0 0 20 4.77 5.07 5.07 0 0 0 19.91 1S18.73.65 16 2.48a13.38 13.38 0 0 0-7 0C6.27.65 5.09 1 5.09 1A5.07 5.07 0 0 0 5 4.77a5.44 5.44 0 0 0-1.5 3.78c0 5.42 3.3 6.61 6.44 7A3.37 3.37 0 0 0 9 18.13V22"
                />
              </svg>
            </a>
            <a :href="`mailto:${SITE.email}`" class="footer-link" :aria-label="'Email'">
              <svg
                class="footer-icon"
                width="16"
                height="16"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                stroke-width="1.75"
                stroke-linecap="round"
                stroke-linejoin="round"
                aria-hidden="true"
              >
                <path
                  d="M4 4h16c1.1 0 2 .9 2 2v12c0 1.1-.9 2-2 2H4c-1.1 0-2-.9-2-2V6c0-1.1.9-2 2-2z"
                />
                <polyline points="22,6 12,13 2,6" />
              </svg>
            </a>
          </div>
        </div>
      </div>
    </div>
  </footer>
</template>

<style scoped>
.site-footer {
  background: var(--app-footer-bg);
  border-top: 1px solid var(--line);
  color: var(--fg-3);
  font-size: var(--text-base);
  padding: var(--sp-8) 0;
  margin-top: var(--sp-12);
}

@media (max-width: 767px) {
  .site-footer .container {
    padding-inline: 20px;
  }
}

.site-footer__inner {
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: var(--sp-6);
}

.footer-copy {
  color: var(--fg-3);
  font-size: var(--text-sm);
  margin: 0;
  line-height: 1.5;
}

.footer-designed {
  display: flex;
  flex-direction: column;
  align-items: center;
  gap: var(--sp-2);
}

.footer-designed__text {
  margin: 0;
  color: var(--fg-3);
  font-size: var(--text-sm);
}

.footer-designed__name {
  font-weight: 700;
  color: var(--fg);
  margin-left: 0.25em;
}

.footer-designed__icons {
  display: flex;
  align-items: center;
  gap: var(--sp-2);
}

.footer-link {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  width: 28px;
  height: 28px;
  border-radius: 50%;
  text-decoration: none;
  transition: background-color 0.14s ease;
}

.footer-link:hover {
  background: var(--tint);
}

.footer-icon {
  display: block;
  color: var(--fg-3);
  opacity: 0.85;
  transition:
    opacity 0.14s ease,
    color 0.14s ease;
}

.footer-link:hover .footer-icon {
  opacity: 1;
  color: var(--primary);
}

@media (max-width: 959.98px) {
  .site-footer__inner {
    flex-direction: column;
    align-items: center;
    text-align: center;
    gap: var(--sp-4);
  }

  .footer-designed {
    flex-direction: row;
    align-items: center;
    gap: var(--sp-2);
  }

  .footer-designed__icons {
    border-left: 1px solid var(--line-strong);
    padding-left: 14px;
  }
}
</style>
